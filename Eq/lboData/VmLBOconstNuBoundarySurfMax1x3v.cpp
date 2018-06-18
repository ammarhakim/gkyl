#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf1x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*2]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[5]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2])+0.6123724356957944*fr[1]*nuVtSq[1]; 

    outr[2] += incr2[2]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*nuVtSq[1]; 

    outl[2] += incr2[2]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*3]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[15]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*(1.369306393762915*fr[12]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[2]*fr[11]+nuVtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[5]); 
    incr2[5] = nuVtSq[1]*(1.369306393762915*fr[12]+0.5477225575051661*fr[11]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[5])+nuVtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[5]); 
    incr2[7] = nuVtSq[0]*(0.6123724356957944*fr[3]-1.060660171779821*fr[7])+0.6123724356957944*nuVtSq[1]*fr[6]; 
    incr2[9] = nuVtSq[0]*(0.6123724356957944*fr[4]-1.060660171779821*fr[9])+0.6123724356957944*nuVtSq[1]*fr[8]; 
    incr2[12] = nuVtSq[0]*((-5.303300858899105*fr[12])+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[2]*fr[11]+nuVtSq[1]*(4.107919181288745*fr[5]-2.371708245126284*fr[1]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[12] += incr2[12]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.369306393762915*fl[12])-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[2]*fl[11]+nuVtSq[1]*((-1.060660171779821*fl[5])-0.6123724356957944*fl[1]); 
    incr2[5] = nuVtSq[1]*((-1.369306393762915*fl[12])-0.5477225575051661*fl[11]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.9486832980505137*fl[5])-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[5])-0.6123724356957944*fl[1]); 
    incr2[7] = nuVtSq[0]*((-1.060660171779821*fl[7])-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[1]*fl[6]; 
    incr2[9] = nuVtSq[0]*((-1.060660171779821*fl[9])-0.6123724356957944*fl[4])-0.6123724356957944*nuVtSq[1]*fl[8]; 
    incr2[12] = nuVtSq[0]*((-5.303300858899105*fl[12])-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[2]*fl[11]+nuVtSq[1]*((-4.107919181288745*fl[5])-2.371708245126284*fl[1]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[12] += incr2[12]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*4]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[35]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*((-1.620185174601965*fr[32])+1.369306393762915*fr[12]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[3]*fr[31]+nuVtSq[1]*(1.369306393762915*fr[20]-1.060660171779821*fr[5]+0.6123724356957944*fr[1])+nuVtSq[2]*(0.6123724356957944*fr[11]-1.060660171779821*fr[19]); 
    incr2[5] = nuVtSq[1]*((-1.620185174601965*fr[32])-0.9486832980505138*fr[19]+1.369306393762915*fr[12]+0.5477225575051661*fr[11]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.537852874200477*fr[31]+1.224744871391589*fr[20]-0.9486832980505137*fr[5]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.369306393762915*fr[20]-1.060660171779821*fr[5]+0.6123724356957944*fr[1])+nuVtSq[3]*(0.537852874200477*fr[11]-0.931588505112178*fr[19]); 
    incr2[7] = nuVtSq[0]*(1.369306393762915*fr[22]-1.060660171779821*fr[7]+0.6123724356957944*fr[3])+0.6123724356957944*nuVtSq[2]*fr[21]+nuVtSq[1]*(0.6123724356957944*fr[6]-1.060660171779821*fr[15]); 
    incr2[9] = nuVtSq[0]*(1.369306393762915*fr[26]-1.060660171779821*fr[9]+0.6123724356957944*fr[4])+0.6123724356957944*nuVtSq[2]*fr[25]+nuVtSq[1]*(0.6123724356957944*fr[8]-1.060660171779821*fr[16]); 
    incr2[12] = nuVtSq[0]*(6.274950199005565*fr[32]-5.303300858899105*fr[12]+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[3]*fr[31]+nuVtSq[1]*((-5.303300858899106*fr[20])+4.107919181288745*fr[5]-2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[19]-2.371708245126284*fr[11]); 
    incr2[15] = nuVtSq[1]*(1.369306393762915*fr[22]+0.5477225575051661*fr[21]-1.060660171779821*fr[7]+0.6123724356957944*fr[3])+0.5378528742004769*nuVtSq[3]*fr[21]+nuVtSq[2]*(0.5477225575051661*fr[6]-0.9486832980505137*fr[15])+nuVtSq[0]*(0.6123724356957944*fr[6]-1.060660171779821*fr[15]); 
    incr2[16] = nuVtSq[1]*(1.369306393762915*fr[26]+0.5477225575051661*fr[25]-1.060660171779821*fr[9]+0.6123724356957944*fr[4])+0.5378528742004769*nuVtSq[3]*fr[25]+nuVtSq[2]*(0.5477225575051661*fr[8]-0.9486832980505137*fr[16])+nuVtSq[0]*(0.6123724356957944*fr[8]-1.060660171779821*fr[16]); 
    incr2[18] = nuVtSq[0]*(0.6123724356957944*fr[10]-1.060660171779821*fr[18])+0.6123724356957944*nuVtSq[1]*fr[17]; 
    incr2[19] = nuVtSq[2]*((-1.620185174601965*fr[32])-0.6776309271789384*fr[19]+1.369306393762915*fr[12]+0.3912303982179757*fr[11]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[1]*(0.5378528742004769*fr[31]+1.224744871391589*fr[20]-0.9486832980505138*fr[5]+0.5477225575051661*fr[1])+nuVtSq[3]*(0.3651483716701107*fr[31]+1.202675588605909*fr[20]-0.931588505112178*fr[5]+0.5378528742004769*fr[1])+nuVtSq[0]*(0.6123724356957944*fr[11]-1.060660171779821*fr[19]); 
    incr2[20] = nuVtSq[1]*(6.274950199005565*fr[32]+3.674234614174766*fr[19]-5.303300858899106*fr[12]-2.121320343559642*fr[11]+4.107919181288746*fr[2]-2.371708245126284*fr[0])+nuVtSq[2]*((-2.08309522448824*fr[31])-4.743416490252569*fr[20]+3.674234614174767*fr[5]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[20])+4.107919181288746*fr[5]-2.371708245126284*fr[1])+nuVtSq[3]*(3.608026765817728*fr[19]-2.08309522448824*fr[11]); 
    incr2[22] = nuVtSq[0]*((-5.303300858899105*fr[22])+4.107919181288745*fr[7]-2.371708245126284*fr[3])-2.371708245126284*nuVtSq[2]*fr[21]+nuVtSq[1]*(4.107919181288745*fr[15]-2.371708245126284*fr[6]); 
    incr2[24] = nuVtSq[0]*(0.6123724356957944*fr[13]-1.060660171779821*fr[24])+0.6123724356957944*nuVtSq[1]*fr[23]; 
    incr2[26] = nuVtSq[0]*((-5.303300858899105*fr[26])+4.107919181288745*fr[9]-2.371708245126284*fr[4])-2.371708245126284*nuVtSq[2]*fr[25]+nuVtSq[1]*(4.107919181288745*fr[16]-2.371708245126284*fr[8]); 
    incr2[29] = nuVtSq[0]*(0.6123724356957944*fr[14]-1.060660171779821*fr[29])+0.6123724356957944*nuVtSq[1]*fr[28]; 
    incr2[32] = nuVtSq[0]*((-14.8492424049175*fr[32])+12.54990039801113*fr[12]-9.721111047611789*fr[2]+5.612486080160912*fr[0])+5.612486080160912*nuVtSq[3]*fr[31]+nuVtSq[1]*(12.54990039801113*fr[20]-9.721111047611789*fr[5]+5.612486080160912*fr[1])+nuVtSq[2]*(5.612486080160912*fr[11]-9.721111047611789*fr[19]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[12] += incr2[12]*rdvSq4r; 
    outr[15] += incr2[15]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[19] += incr2[19]*rdvSq4r; 
    outr[20] += incr2[20]*rdvSq4r; 
    outr[22] += incr2[22]*rdvSq4r; 
    outr[24] += incr2[24]*rdvSq4r; 
    outr[26] += incr2[26]*rdvSq4r; 
    outr[29] += incr2[29]*rdvSq4r; 
    outr[32] += incr2[32]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.620185174601965*fl[32])-1.369306393762915*fl[12]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[3]*fl[31]+nuVtSq[1]*((-1.369306393762915*fl[20])-1.060660171779821*fl[5]-0.6123724356957944*fl[1])+nuVtSq[2]*((-1.060660171779821*fl[19])-0.6123724356957944*fl[11]); 
    incr2[5] = nuVtSq[1]*((-1.620185174601965*fl[32])-0.9486832980505138*fl[19]-1.369306393762915*fl[12]-0.5477225575051661*fl[11]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.537852874200477*fl[31])-1.224744871391589*fl[20]-0.9486832980505137*fl[5]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[20])-1.060660171779821*fl[5]-0.6123724356957944*fl[1])+nuVtSq[3]*((-0.931588505112178*fl[19])-0.537852874200477*fl[11]); 
    incr2[7] = nuVtSq[0]*((-1.369306393762915*fl[22])-1.060660171779821*fl[7]-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[2]*fl[21]+nuVtSq[1]*((-1.060660171779821*fl[15])-0.6123724356957944*fl[6]); 
    incr2[9] = nuVtSq[0]*((-1.369306393762915*fl[26])-1.060660171779821*fl[9]-0.6123724356957944*fl[4])-0.6123724356957944*nuVtSq[2]*fl[25]+nuVtSq[1]*((-1.060660171779821*fl[16])-0.6123724356957944*fl[8]); 
    incr2[12] = nuVtSq[0]*((-6.274950199005565*fl[32])-5.303300858899105*fl[12]-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[3]*fl[31]+nuVtSq[1]*((-5.303300858899106*fl[20])-4.107919181288745*fl[5]-2.371708245126284*fl[1])+nuVtSq[2]*((-4.107919181288745*fl[19])-2.371708245126284*fl[11]); 
    incr2[15] = nuVtSq[1]*((-1.369306393762915*fl[22])-0.5477225575051661*fl[21]-1.060660171779821*fl[7]-0.6123724356957944*fl[3])-0.5378528742004769*nuVtSq[3]*fl[21]+nuVtSq[2]*((-0.9486832980505137*fl[15])-0.5477225575051661*fl[6])+nuVtSq[0]*((-1.060660171779821*fl[15])-0.6123724356957944*fl[6]); 
    incr2[16] = nuVtSq[1]*((-1.369306393762915*fl[26])-0.5477225575051661*fl[25]-1.060660171779821*fl[9]-0.6123724356957944*fl[4])-0.5378528742004769*nuVtSq[3]*fl[25]+nuVtSq[2]*((-0.9486832980505137*fl[16])-0.5477225575051661*fl[8])+nuVtSq[0]*((-1.060660171779821*fl[16])-0.6123724356957944*fl[8]); 
    incr2[18] = nuVtSq[0]*((-1.060660171779821*fl[18])-0.6123724356957944*fl[10])-0.6123724356957944*nuVtSq[1]*fl[17]; 
    incr2[19] = nuVtSq[2]*((-1.620185174601965*fl[32])-0.6776309271789384*fl[19]-1.369306393762915*fl[12]-0.3912303982179757*fl[11]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[3]*((-0.3651483716701107*fl[31])-1.202675588605909*fl[20]-0.931588505112178*fl[5]-0.5378528742004769*fl[1])+nuVtSq[1]*((-0.5378528742004769*fl[31])-1.224744871391589*fl[20]-0.9486832980505138*fl[5]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[19])-0.6123724356957944*fl[11]); 
    incr2[20] = nuVtSq[1]*((-6.274950199005565*fl[32])-3.674234614174766*fl[19]-5.303300858899106*fl[12]-2.121320343559642*fl[11]-4.107919181288746*fl[2]-2.371708245126284*fl[0])+nuVtSq[2]*((-2.08309522448824*fl[31])-4.743416490252569*fl[20]-3.674234614174767*fl[5]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[20])-4.107919181288746*fl[5]-2.371708245126284*fl[1])+nuVtSq[3]*((-3.608026765817728*fl[19])-2.08309522448824*fl[11]); 
    incr2[22] = nuVtSq[0]*((-5.303300858899105*fl[22])-4.107919181288745*fl[7]-2.371708245126284*fl[3])-2.371708245126284*nuVtSq[2]*fl[21]+nuVtSq[1]*((-4.107919181288745*fl[15])-2.371708245126284*fl[6]); 
    incr2[24] = nuVtSq[0]*((-1.060660171779821*fl[24])-0.6123724356957944*fl[13])-0.6123724356957944*nuVtSq[1]*fl[23]; 
    incr2[26] = nuVtSq[0]*((-5.303300858899105*fl[26])-4.107919181288745*fl[9]-2.371708245126284*fl[4])-2.371708245126284*nuVtSq[2]*fl[25]+nuVtSq[1]*((-4.107919181288745*fl[16])-2.371708245126284*fl[8]); 
    incr2[29] = nuVtSq[0]*((-1.060660171779821*fl[29])-0.6123724356957944*fl[14])-0.6123724356957944*nuVtSq[1]*fl[28]; 
    incr2[32] = nuVtSq[0]*((-14.8492424049175*fl[32])-12.54990039801113*fl[12]-9.721111047611789*fl[2]-5.612486080160912*fl[0])-5.612486080160912*nuVtSq[3]*fl[31]+nuVtSq[1]*((-12.54990039801113*fl[20])-9.721111047611789*fl[5]-5.612486080160912*fl[1])+nuVtSq[2]*((-9.721111047611789*fl[19])-5.612486080160912*fl[11]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[12] += incr2[12]*rdvSq4l; 
    outl[15] += incr2[15]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[19] += incr2[19]*rdvSq4l; 
    outl[20] += incr2[20]*rdvSq4l; 
    outl[22] += incr2[22]*rdvSq4l; 
    outl[24] += incr2[24]*rdvSq4l; 
    outl[26] += incr2[26]*rdvSq4l; 
    outl[29] += incr2[29]*rdvSq4l; 
    outl[32] += incr2[32]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*2]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[5]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[3])+0.6123724356957944*fr[1]*nuVtSq[1]; 

    outr[3] += incr2[3]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*nuVtSq[1]; 

    outl[3] += incr2[3]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*3]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[15]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*(1.369306393762915*fr[13]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[2]*fr[11]+nuVtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[6]); 
    incr2[6] = nuVtSq[1]*(1.369306393762915*fr[13]+0.5477225575051661*fr[11]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[6])+nuVtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[6]); 
    incr2[7] = nuVtSq[0]*(0.6123724356957944*fr[2]-1.060660171779821*fr[7])+0.6123724356957944*nuVtSq[1]*fr[5]; 
    incr2[10] = nuVtSq[0]*(0.6123724356957944*fr[4]-1.060660171779821*fr[10])+0.6123724356957944*nuVtSq[1]*fr[8]; 
    incr2[13] = nuVtSq[0]*((-5.303300858899105*fr[13])+4.107919181288745*fr[3]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[2]*fr[11]+nuVtSq[1]*(4.107919181288745*fr[6]-2.371708245126284*fr[1]); 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[13] += incr2[13]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.369306393762915*fl[13])-1.060660171779821*fl[3]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[2]*fl[11]+nuVtSq[1]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[1]); 
    incr2[6] = nuVtSq[1]*((-1.369306393762915*fl[13])-0.5477225575051661*fl[11]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.9486832980505137*fl[6])-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[1]); 
    incr2[7] = nuVtSq[0]*((-1.060660171779821*fl[7])-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[1]*fl[5]; 
    incr2[10] = nuVtSq[0]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[4])-0.6123724356957944*nuVtSq[1]*fl[8]; 
    incr2[13] = nuVtSq[0]*((-5.303300858899105*fl[13])-4.107919181288745*fl[3]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[2]*fl[11]+nuVtSq[1]*((-4.107919181288745*fl[6])-2.371708245126284*fl[1]); 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[13] += incr2[13]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*4]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[35]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*((-1.620185174601965*fr[33])+1.369306393762915*fr[13]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[3]*fr[31]+nuVtSq[1]*(1.369306393762915*fr[23]-1.060660171779821*fr[6]+0.6123724356957944*fr[1])+nuVtSq[2]*(0.6123724356957944*fr[11]-1.060660171779821*fr[21]); 
    incr2[6] = nuVtSq[1]*((-1.620185174601965*fr[33])-0.9486832980505138*fr[21]+1.369306393762915*fr[13]+0.5477225575051661*fr[11]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.537852874200477*fr[31]+1.224744871391589*fr[23]-0.9486832980505137*fr[6]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.369306393762915*fr[23]-1.060660171779821*fr[6]+0.6123724356957944*fr[1])+nuVtSq[3]*(0.537852874200477*fr[11]-0.931588505112178*fr[21]); 
    incr2[7] = nuVtSq[0]*(1.369306393762915*fr[24]-1.060660171779821*fr[7]+0.6123724356957944*fr[2])+0.6123724356957944*nuVtSq[2]*fr[19]+nuVtSq[1]*(0.6123724356957944*fr[5]-1.060660171779821*fr[15]); 
    incr2[10] = nuVtSq[0]*(1.369306393762915*fr[27]-1.060660171779821*fr[10]+0.6123724356957944*fr[4])+0.6123724356957944*nuVtSq[2]*fr[25]+nuVtSq[1]*(0.6123724356957944*fr[8]-1.060660171779821*fr[17]); 
    incr2[13] = nuVtSq[0]*(6.274950199005565*fr[33]-5.303300858899105*fr[13]+4.107919181288745*fr[3]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[3]*fr[31]+nuVtSq[1]*((-5.303300858899106*fr[23])+4.107919181288745*fr[6]-2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[21]-2.371708245126284*fr[11]); 
    incr2[15] = nuVtSq[1]*(1.369306393762915*fr[24]+0.5477225575051661*fr[19]-1.060660171779821*fr[7]+0.6123724356957944*fr[2])+0.5378528742004769*nuVtSq[3]*fr[19]+nuVtSq[2]*(0.5477225575051661*fr[5]-0.9486832980505137*fr[15])+nuVtSq[0]*(0.6123724356957944*fr[5]-1.060660171779821*fr[15]); 
    incr2[17] = nuVtSq[1]*(1.369306393762915*fr[27]+0.5477225575051661*fr[25]-1.060660171779821*fr[10]+0.6123724356957944*fr[4])+0.5378528742004769*nuVtSq[3]*fr[25]+nuVtSq[2]*(0.5477225575051661*fr[8]-0.9486832980505137*fr[17])+nuVtSq[0]*(0.6123724356957944*fr[8]-1.060660171779821*fr[17]); 
    incr2[18] = nuVtSq[0]*(0.6123724356957944*fr[9]-1.060660171779821*fr[18])+0.6123724356957944*nuVtSq[1]*fr[16]; 
    incr2[21] = nuVtSq[2]*((-1.620185174601965*fr[33])-0.6776309271789384*fr[21]+1.369306393762915*fr[13]+0.3912303982179757*fr[11]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[1]*(0.5378528742004769*fr[31]+1.224744871391589*fr[23]-0.9486832980505138*fr[6]+0.5477225575051661*fr[1])+nuVtSq[3]*(0.3651483716701107*fr[31]+1.202675588605909*fr[23]-0.931588505112178*fr[6]+0.5378528742004769*fr[1])+nuVtSq[0]*(0.6123724356957944*fr[11]-1.060660171779821*fr[21]); 
    incr2[22] = nuVtSq[0]*(0.6123724356957944*fr[12]-1.060660171779821*fr[22])+0.6123724356957944*nuVtSq[1]*fr[20]; 
    incr2[23] = nuVtSq[1]*(6.274950199005565*fr[33]+3.674234614174766*fr[21]-5.303300858899106*fr[13]-2.121320343559642*fr[11]+4.107919181288746*fr[3]-2.371708245126284*fr[0])+nuVtSq[2]*((-2.08309522448824*fr[31])-4.743416490252569*fr[23]+3.674234614174767*fr[6]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[23])+4.107919181288746*fr[6]-2.371708245126284*fr[1])+nuVtSq[3]*(3.608026765817728*fr[21]-2.08309522448824*fr[11]); 
    incr2[24] = nuVtSq[0]*((-5.303300858899105*fr[24])+4.107919181288745*fr[7]-2.371708245126284*fr[2])-2.371708245126284*nuVtSq[2]*fr[19]+nuVtSq[1]*(4.107919181288745*fr[15]-2.371708245126284*fr[5]); 
    incr2[27] = nuVtSq[0]*((-5.303300858899105*fr[27])+4.107919181288745*fr[10]-2.371708245126284*fr[4])-2.371708245126284*nuVtSq[2]*fr[25]+nuVtSq[1]*(4.107919181288745*fr[17]-2.371708245126284*fr[8]); 
    incr2[30] = nuVtSq[0]*(0.6123724356957944*fr[14]-1.060660171779821*fr[30])+0.6123724356957944*nuVtSq[1]*fr[28]; 
    incr2[33] = nuVtSq[0]*((-14.8492424049175*fr[33])+12.54990039801113*fr[13]-9.721111047611789*fr[3]+5.612486080160912*fr[0])+5.612486080160912*nuVtSq[3]*fr[31]+nuVtSq[1]*(12.54990039801113*fr[23]-9.721111047611789*fr[6]+5.612486080160912*fr[1])+nuVtSq[2]*(5.612486080160912*fr[11]-9.721111047611789*fr[21]); 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[13] += incr2[13]*rdvSq4r; 
    outr[15] += incr2[15]*rdvSq4r; 
    outr[17] += incr2[17]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[21] += incr2[21]*rdvSq4r; 
    outr[22] += incr2[22]*rdvSq4r; 
    outr[23] += incr2[23]*rdvSq4r; 
    outr[24] += incr2[24]*rdvSq4r; 
    outr[27] += incr2[27]*rdvSq4r; 
    outr[30] += incr2[30]*rdvSq4r; 
    outr[33] += incr2[33]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.620185174601965*fl[33])-1.369306393762915*fl[13]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[3]*fl[31]+nuVtSq[1]*((-1.369306393762915*fl[23])-1.060660171779821*fl[6]-0.6123724356957944*fl[1])+nuVtSq[2]*((-1.060660171779821*fl[21])-0.6123724356957944*fl[11]); 
    incr2[6] = nuVtSq[1]*((-1.620185174601965*fl[33])-0.9486832980505138*fl[21]-1.369306393762915*fl[13]-0.5477225575051661*fl[11]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.537852874200477*fl[31])-1.224744871391589*fl[23]-0.9486832980505137*fl[6]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[23])-1.060660171779821*fl[6]-0.6123724356957944*fl[1])+nuVtSq[3]*((-0.931588505112178*fl[21])-0.537852874200477*fl[11]); 
    incr2[7] = nuVtSq[0]*((-1.369306393762915*fl[24])-1.060660171779821*fl[7]-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[2]*fl[19]+nuVtSq[1]*((-1.060660171779821*fl[15])-0.6123724356957944*fl[5]); 
    incr2[10] = nuVtSq[0]*((-1.369306393762915*fl[27])-1.060660171779821*fl[10]-0.6123724356957944*fl[4])-0.6123724356957944*nuVtSq[2]*fl[25]+nuVtSq[1]*((-1.060660171779821*fl[17])-0.6123724356957944*fl[8]); 
    incr2[13] = nuVtSq[0]*((-6.274950199005565*fl[33])-5.303300858899105*fl[13]-4.107919181288745*fl[3]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[3]*fl[31]+nuVtSq[1]*((-5.303300858899106*fl[23])-4.107919181288745*fl[6]-2.371708245126284*fl[1])+nuVtSq[2]*((-4.107919181288745*fl[21])-2.371708245126284*fl[11]); 
    incr2[15] = nuVtSq[1]*((-1.369306393762915*fl[24])-0.5477225575051661*fl[19]-1.060660171779821*fl[7]-0.6123724356957944*fl[2])-0.5378528742004769*nuVtSq[3]*fl[19]+nuVtSq[2]*((-0.9486832980505137*fl[15])-0.5477225575051661*fl[5])+nuVtSq[0]*((-1.060660171779821*fl[15])-0.6123724356957944*fl[5]); 
    incr2[17] = nuVtSq[1]*((-1.369306393762915*fl[27])-0.5477225575051661*fl[25]-1.060660171779821*fl[10]-0.6123724356957944*fl[4])-0.5378528742004769*nuVtSq[3]*fl[25]+nuVtSq[2]*((-0.9486832980505137*fl[17])-0.5477225575051661*fl[8])+nuVtSq[0]*((-1.060660171779821*fl[17])-0.6123724356957944*fl[8]); 
    incr2[18] = nuVtSq[0]*((-1.060660171779821*fl[18])-0.6123724356957944*fl[9])-0.6123724356957944*nuVtSq[1]*fl[16]; 
    incr2[21] = nuVtSq[2]*((-1.620185174601965*fl[33])-0.6776309271789384*fl[21]-1.369306393762915*fl[13]-0.3912303982179757*fl[11]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[3]*((-0.3651483716701107*fl[31])-1.202675588605909*fl[23]-0.931588505112178*fl[6]-0.5378528742004769*fl[1])+nuVtSq[1]*((-0.5378528742004769*fl[31])-1.224744871391589*fl[23]-0.9486832980505138*fl[6]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[21])-0.6123724356957944*fl[11]); 
    incr2[22] = nuVtSq[0]*((-1.060660171779821*fl[22])-0.6123724356957944*fl[12])-0.6123724356957944*nuVtSq[1]*fl[20]; 
    incr2[23] = nuVtSq[1]*((-6.274950199005565*fl[33])-3.674234614174766*fl[21]-5.303300858899106*fl[13]-2.121320343559642*fl[11]-4.107919181288746*fl[3]-2.371708245126284*fl[0])+nuVtSq[2]*((-2.08309522448824*fl[31])-4.743416490252569*fl[23]-3.674234614174767*fl[6]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[23])-4.107919181288746*fl[6]-2.371708245126284*fl[1])+nuVtSq[3]*((-3.608026765817728*fl[21])-2.08309522448824*fl[11]); 
    incr2[24] = nuVtSq[0]*((-5.303300858899105*fl[24])-4.107919181288745*fl[7]-2.371708245126284*fl[2])-2.371708245126284*nuVtSq[2]*fl[19]+nuVtSq[1]*((-4.107919181288745*fl[15])-2.371708245126284*fl[5]); 
    incr2[27] = nuVtSq[0]*((-5.303300858899105*fl[27])-4.107919181288745*fl[10]-2.371708245126284*fl[4])-2.371708245126284*nuVtSq[2]*fl[25]+nuVtSq[1]*((-4.107919181288745*fl[17])-2.371708245126284*fl[8]); 
    incr2[30] = nuVtSq[0]*((-1.060660171779821*fl[30])-0.6123724356957944*fl[14])-0.6123724356957944*nuVtSq[1]*fl[28]; 
    incr2[33] = nuVtSq[0]*((-14.8492424049175*fl[33])-12.54990039801113*fl[13]-9.721111047611789*fl[3]-5.612486080160912*fl[0])-5.612486080160912*nuVtSq[3]*fl[31]+nuVtSq[1]*((-12.54990039801113*fl[23])-9.721111047611789*fl[6]-5.612486080160912*fl[1])+nuVtSq[2]*((-9.721111047611789*fl[21])-5.612486080160912*fl[11]); 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[13] += incr2[13]*rdvSq4l; 
    outl[15] += incr2[15]*rdvSq4l; 
    outl[17] += incr2[17]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[21] += incr2[21]*rdvSq4l; 
    outl[22] += incr2[22]*rdvSq4l; 
    outl[23] += incr2[23]*rdvSq4l; 
    outl[24] += incr2[24]*rdvSq4l; 
    outl[27] += incr2[27]*rdvSq4l; 
    outl[30] += incr2[30]*rdvSq4l; 
    outl[33] += incr2[33]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*2]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4r = 4.0/(dxvr[3]*dxvr[3]); 

  double incr2[5]; 
  if (idxr[3] == 1) {

    incr2[4] = nuVtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[4])+0.6123724356957944*fr[1]*nuVtSq[1]; 

    outr[4] += incr2[4]*rdvSq4r; 

  } else {

    incr2[4] = nuVtSq[0]*((-1.060660171779821*fl[4])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*nuVtSq[1]; 

    outl[4] += incr2[4]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*3]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4r = 4.0/(dxvr[3]*dxvr[3]); 

  double incr2[15]; 
  if (idxr[3] == 1) {

    incr2[4] = nuVtSq[0]*(1.369306393762915*fr[14]-1.060660171779821*fr[4]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[2]*fr[11]+nuVtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[8]); 
    incr2[8] = nuVtSq[1]*(1.369306393762915*fr[14]+0.5477225575051661*fr[11]-1.060660171779821*fr[4]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[8])+nuVtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[8]); 
    incr2[9] = nuVtSq[0]*(0.6123724356957944*fr[2]-1.060660171779821*fr[9])+0.6123724356957944*nuVtSq[1]*fr[5]; 
    incr2[10] = nuVtSq[0]*(0.6123724356957944*fr[3]-1.060660171779821*fr[10])+0.6123724356957944*nuVtSq[1]*fr[6]; 
    incr2[14] = nuVtSq[0]*((-5.303300858899105*fr[14])+4.107919181288745*fr[4]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[2]*fr[11]+nuVtSq[1]*(4.107919181288745*fr[8]-2.371708245126284*fr[1]); 

    outr[4] += incr2[4]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 

  } else {

    incr2[4] = nuVtSq[0]*((-1.369306393762915*fl[14])-1.060660171779821*fl[4]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[2]*fl[11]+nuVtSq[1]*((-1.060660171779821*fl[8])-0.6123724356957944*fl[1]); 
    incr2[8] = nuVtSq[1]*((-1.369306393762915*fl[14])-0.5477225575051661*fl[11]-1.060660171779821*fl[4]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.9486832980505137*fl[8])-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[8])-0.6123724356957944*fl[1]); 
    incr2[9] = nuVtSq[0]*((-1.060660171779821*fl[9])-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[1]*fl[5]; 
    incr2[10] = nuVtSq[0]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[1]*fl[6]; 
    incr2[14] = nuVtSq[0]*((-5.303300858899105*fl[14])-4.107919181288745*fl[4]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[2]*fl[11]+nuVtSq[1]*((-4.107919181288745*fl[8])-2.371708245126284*fl[1]); 

    outl[4] += incr2[4]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vMax_VZ_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[3*4]: bulk velocity (in 3 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4r = 4.0/(dxvr[3]*dxvr[3]); 

  double incr2[35]; 
  if (idxr[3] == 1) {

    incr2[4] = nuVtSq[0]*((-1.620185174601965*fr[34])+1.369306393762915*fr[14]-1.060660171779821*fr[4]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[3]*fr[31]+nuVtSq[1]*(1.369306393762915*fr[28]-1.060660171779821*fr[8]+0.6123724356957944*fr[1])+nuVtSq[2]*(0.6123724356957944*fr[11]-1.060660171779821*fr[25]); 
    incr2[8] = nuVtSq[1]*((-1.620185174601965*fr[34])-0.9486832980505138*fr[25]+1.369306393762915*fr[14]+0.5477225575051661*fr[11]-1.060660171779821*fr[4]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.537852874200477*fr[31]+1.224744871391589*fr[28]-0.9486832980505137*fr[8]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.369306393762915*fr[28]-1.060660171779821*fr[8]+0.6123724356957944*fr[1])+nuVtSq[3]*(0.537852874200477*fr[11]-0.931588505112178*fr[25]); 
    incr2[9] = nuVtSq[0]*(1.369306393762915*fr[29]-1.060660171779821*fr[9]+0.6123724356957944*fr[2])+0.6123724356957944*nuVtSq[2]*fr[19]+nuVtSq[1]*(0.6123724356957944*fr[5]-1.060660171779821*fr[16]); 
    incr2[10] = nuVtSq[0]*(1.369306393762915*fr[30]-1.060660171779821*fr[10]+0.6123724356957944*fr[3])+0.6123724356957944*nuVtSq[2]*fr[21]+nuVtSq[1]*(0.6123724356957944*fr[6]-1.060660171779821*fr[17]); 
    incr2[14] = nuVtSq[0]*(6.274950199005565*fr[34]-5.303300858899105*fr[14]+4.107919181288745*fr[4]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[3]*fr[31]+nuVtSq[1]*((-5.303300858899106*fr[28])+4.107919181288745*fr[8]-2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[25]-2.371708245126284*fr[11]); 
    incr2[16] = nuVtSq[1]*(1.369306393762915*fr[29]+0.5477225575051661*fr[19]-1.060660171779821*fr[9]+0.6123724356957944*fr[2])+0.5378528742004769*nuVtSq[3]*fr[19]+nuVtSq[2]*(0.5477225575051661*fr[5]-0.9486832980505137*fr[16])+nuVtSq[0]*(0.6123724356957944*fr[5]-1.060660171779821*fr[16]); 
    incr2[17] = nuVtSq[1]*(1.369306393762915*fr[30]+0.5477225575051661*fr[21]-1.060660171779821*fr[10]+0.6123724356957944*fr[3])+0.5378528742004769*nuVtSq[3]*fr[21]+nuVtSq[2]*(0.5477225575051661*fr[6]-0.9486832980505137*fr[17])+nuVtSq[0]*(0.6123724356957944*fr[6]-1.060660171779821*fr[17]); 
    incr2[18] = nuVtSq[0]*(0.6123724356957944*fr[7]-1.060660171779821*fr[18])+0.6123724356957944*nuVtSq[1]*fr[15]; 
    incr2[25] = nuVtSq[2]*((-1.620185174601965*fr[34])-0.6776309271789384*fr[25]+1.369306393762915*fr[14]+0.3912303982179757*fr[11]-1.060660171779821*fr[4]+0.6123724356957944*fr[0])+nuVtSq[1]*(0.5378528742004769*fr[31]+1.224744871391589*fr[28]-0.9486832980505138*fr[8]+0.5477225575051661*fr[1])+nuVtSq[3]*(0.3651483716701107*fr[31]+1.202675588605909*fr[28]-0.931588505112178*fr[8]+0.5378528742004769*fr[1])+nuVtSq[0]*(0.6123724356957944*fr[11]-1.060660171779821*fr[25]); 
    incr2[26] = nuVtSq[0]*(0.6123724356957944*fr[12]-1.060660171779821*fr[26])+0.6123724356957944*nuVtSq[1]*fr[20]; 
    incr2[27] = nuVtSq[0]*(0.6123724356957944*fr[13]-1.060660171779821*fr[27])+0.6123724356957944*nuVtSq[1]*fr[23]; 
    incr2[28] = nuVtSq[1]*(6.274950199005565*fr[34]+3.674234614174766*fr[25]-5.303300858899106*fr[14]-2.121320343559642*fr[11]+4.107919181288746*fr[4]-2.371708245126284*fr[0])+nuVtSq[2]*((-2.08309522448824*fr[31])-4.743416490252569*fr[28]+3.674234614174767*fr[8]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[28])+4.107919181288746*fr[8]-2.371708245126284*fr[1])+nuVtSq[3]*(3.608026765817728*fr[25]-2.08309522448824*fr[11]); 
    incr2[29] = nuVtSq[0]*((-5.303300858899105*fr[29])+4.107919181288745*fr[9]-2.371708245126284*fr[2])-2.371708245126284*nuVtSq[2]*fr[19]+nuVtSq[1]*(4.107919181288745*fr[16]-2.371708245126284*fr[5]); 
    incr2[30] = nuVtSq[0]*((-5.303300858899105*fr[30])+4.107919181288745*fr[10]-2.371708245126284*fr[3])-2.371708245126284*nuVtSq[2]*fr[21]+nuVtSq[1]*(4.107919181288745*fr[17]-2.371708245126284*fr[6]); 
    incr2[34] = nuVtSq[0]*((-14.8492424049175*fr[34])+12.54990039801113*fr[14]-9.721111047611789*fr[4]+5.612486080160912*fr[0])+5.612486080160912*nuVtSq[3]*fr[31]+nuVtSq[1]*(12.54990039801113*fr[28]-9.721111047611789*fr[8]+5.612486080160912*fr[1])+nuVtSq[2]*(5.612486080160912*fr[11]-9.721111047611789*fr[25]); 

    outr[4] += incr2[4]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[17] += incr2[17]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[25] += incr2[25]*rdvSq4r; 
    outr[26] += incr2[26]*rdvSq4r; 
    outr[27] += incr2[27]*rdvSq4r; 
    outr[28] += incr2[28]*rdvSq4r; 
    outr[29] += incr2[29]*rdvSq4r; 
    outr[30] += incr2[30]*rdvSq4r; 
    outr[34] += incr2[34]*rdvSq4r; 

  } else {

    incr2[4] = nuVtSq[0]*((-1.620185174601965*fl[34])-1.369306393762915*fl[14]-1.060660171779821*fl[4]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[3]*fl[31]+nuVtSq[1]*((-1.369306393762915*fl[28])-1.060660171779821*fl[8]-0.6123724356957944*fl[1])+nuVtSq[2]*((-1.060660171779821*fl[25])-0.6123724356957944*fl[11]); 
    incr2[8] = nuVtSq[1]*((-1.620185174601965*fl[34])-0.9486832980505138*fl[25]-1.369306393762915*fl[14]-0.5477225575051661*fl[11]-1.060660171779821*fl[4]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.537852874200477*fl[31])-1.224744871391589*fl[28]-0.9486832980505137*fl[8]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[28])-1.060660171779821*fl[8]-0.6123724356957944*fl[1])+nuVtSq[3]*((-0.931588505112178*fl[25])-0.537852874200477*fl[11]); 
    incr2[9] = nuVtSq[0]*((-1.369306393762915*fl[29])-1.060660171779821*fl[9]-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[2]*fl[19]+nuVtSq[1]*((-1.060660171779821*fl[16])-0.6123724356957944*fl[5]); 
    incr2[10] = nuVtSq[0]*((-1.369306393762915*fl[30])-1.060660171779821*fl[10]-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[2]*fl[21]+nuVtSq[1]*((-1.060660171779821*fl[17])-0.6123724356957944*fl[6]); 
    incr2[14] = nuVtSq[0]*((-6.274950199005565*fl[34])-5.303300858899105*fl[14]-4.107919181288745*fl[4]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[3]*fl[31]+nuVtSq[1]*((-5.303300858899106*fl[28])-4.107919181288745*fl[8]-2.371708245126284*fl[1])+nuVtSq[2]*((-4.107919181288745*fl[25])-2.371708245126284*fl[11]); 
    incr2[16] = nuVtSq[1]*((-1.369306393762915*fl[29])-0.5477225575051661*fl[19]-1.060660171779821*fl[9]-0.6123724356957944*fl[2])-0.5378528742004769*nuVtSq[3]*fl[19]+nuVtSq[2]*((-0.9486832980505137*fl[16])-0.5477225575051661*fl[5])+nuVtSq[0]*((-1.060660171779821*fl[16])-0.6123724356957944*fl[5]); 
    incr2[17] = nuVtSq[1]*((-1.369306393762915*fl[30])-0.5477225575051661*fl[21]-1.060660171779821*fl[10]-0.6123724356957944*fl[3])-0.5378528742004769*nuVtSq[3]*fl[21]+nuVtSq[2]*((-0.9486832980505137*fl[17])-0.5477225575051661*fl[6])+nuVtSq[0]*((-1.060660171779821*fl[17])-0.6123724356957944*fl[6]); 
    incr2[18] = nuVtSq[0]*((-1.060660171779821*fl[18])-0.6123724356957944*fl[7])-0.6123724356957944*nuVtSq[1]*fl[15]; 
    incr2[25] = nuVtSq[2]*((-1.620185174601965*fl[34])-0.6776309271789384*fl[25]-1.369306393762915*fl[14]-0.3912303982179757*fl[11]-1.060660171779821*fl[4]-0.6123724356957944*fl[0])+nuVtSq[3]*((-0.3651483716701107*fl[31])-1.202675588605909*fl[28]-0.931588505112178*fl[8]-0.5378528742004769*fl[1])+nuVtSq[1]*((-0.5378528742004769*fl[31])-1.224744871391589*fl[28]-0.9486832980505138*fl[8]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[25])-0.6123724356957944*fl[11]); 
    incr2[26] = nuVtSq[0]*((-1.060660171779821*fl[26])-0.6123724356957944*fl[12])-0.6123724356957944*nuVtSq[1]*fl[20]; 
    incr2[27] = nuVtSq[0]*((-1.060660171779821*fl[27])-0.6123724356957944*fl[13])-0.6123724356957944*nuVtSq[1]*fl[23]; 
    incr2[28] = nuVtSq[1]*((-6.274950199005565*fl[34])-3.674234614174766*fl[25]-5.303300858899106*fl[14]-2.121320343559642*fl[11]-4.107919181288746*fl[4]-2.371708245126284*fl[0])+nuVtSq[2]*((-2.08309522448824*fl[31])-4.743416490252569*fl[28]-3.674234614174767*fl[8]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[28])-4.107919181288746*fl[8]-2.371708245126284*fl[1])+nuVtSq[3]*((-3.608026765817728*fl[25])-2.08309522448824*fl[11]); 
    incr2[29] = nuVtSq[0]*((-5.303300858899105*fl[29])-4.107919181288745*fl[9]-2.371708245126284*fl[2])-2.371708245126284*nuVtSq[2]*fl[19]+nuVtSq[1]*((-4.107919181288745*fl[16])-2.371708245126284*fl[5]); 
    incr2[30] = nuVtSq[0]*((-5.303300858899105*fl[30])-4.107919181288745*fl[10]-2.371708245126284*fl[3])-2.371708245126284*nuVtSq[2]*fl[21]+nuVtSq[1]*((-4.107919181288745*fl[17])-2.371708245126284*fl[6]); 
    incr2[34] = nuVtSq[0]*((-14.8492424049175*fl[34])-12.54990039801113*fl[14]-9.721111047611789*fl[4]-5.612486080160912*fl[0])-5.612486080160912*nuVtSq[3]*fl[31]+nuVtSq[1]*((-12.54990039801113*fl[28])-9.721111047611789*fl[8]-5.612486080160912*fl[1])+nuVtSq[2]*((-9.721111047611789*fl[25])-5.612486080160912*fl[11]); 

    outl[4] += incr2[4]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[17] += incr2[17]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[25] += incr2[25]*rdvSq4l; 
    outl[26] += incr2[26]*rdvSq4l; 
    outl[27] += incr2[27]*rdvSq4l; 
    outl[28] += incr2[28]*rdvSq4l; 
    outl[29] += incr2[29]*rdvSq4l; 
    outl[30] += incr2[30]*rdvSq4l; 
    outl[34] += incr2[34]*rdvSq4l; 

  }
  return 0.0; 
} 
