#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[8]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[1]*(1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[0]*(1.060660171779821*fr[2]+0.6123724356957944*fr[0]); 
    incr2[4] = nuVtSq[0]*(1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[1]*(1.060660171779821*fr[2]+0.6123724356957944*fr[0]); 
    incr2[6] = nuVtSq[1]*(1.060660171779821*fr[7]+0.6123724356957944*fr[5])+nuVtSq[0]*(1.060660171779821*fr[6]+0.6123724356957944*fr[3]); 
    incr2[7] = nuVtSq[0]*(1.060660171779821*fr[7]+0.6123724356957944*fr[5])+nuVtSq[1]*(1.060660171779821*fr[6]+0.6123724356957944*fr[3]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[7] += incr2[7]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[1]*(1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[0]*(1.060660171779821*fl[2]-0.6123724356957944*fl[0]); 
    incr2[4] = nuVtSq[0]*(1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[1]*(1.060660171779821*fl[2]-0.6123724356957944*fl[0]); 
    incr2[6] = nuVtSq[1]*(1.060660171779821*fl[7]-0.6123724356957944*fl[5])+nuVtSq[0]*(1.060660171779821*fl[6]-0.6123724356957944*fl[3]); 
    incr2[7] = nuVtSq[0]*(1.060660171779821*fl[7]-0.6123724356957944*fl[5])+nuVtSq[1]*(1.060660171779821*fl[6]-0.6123724356957944*fl[3]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[7] += incr2[7]*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[20]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[1]*(1.369306393762915*fr[12]+1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[2]*(1.060660171779821*fr[11]+0.6123724356957944*fr[7])+nuVtSq[0]*(1.369306393762915*fr[8]+1.060660171779821*fr[2]+0.6123724356957944*fr[0]); 
    incr2[4] = nuVtSq[0]*(1.369306393762915*fr[12]+1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[2]*(1.224744871391589*fr[12]+0.9486832980505137*fr[4]+0.5477225575051661*fr[1])+nuVtSq[1]*(0.9486832980505138*fr[11]+1.369306393762915*fr[8]+0.5477225575051661*fr[7]+1.060660171779821*fr[2]+0.6123724356957944*fr[0]); 
    incr2[6] = nuVtSq[1]*(1.369306393762915*fr[18]+1.060660171779821*fr[10]+0.6123724356957944*fr[5])+nuVtSq[2]*(1.060660171779821*fr[17]+0.6123724356957944*fr[13])+nuVtSq[0]*(1.369306393762915*fr[14]+1.060660171779821*fr[6]+0.6123724356957944*fr[3]); 
    incr2[8] = nuVtSq[1]*(5.303300858899106*fr[12]+4.107919181288745*fr[4]+2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[11]+2.371708245126284*fr[7])+nuVtSq[0]*(5.303300858899105*fr[8]+4.107919181288745*fr[2]+2.371708245126284*fr[0]); 
    incr2[10] = nuVtSq[0]*(1.369306393762915*fr[18]+1.060660171779821*fr[10]+0.6123724356957944*fr[5])+nuVtSq[2]*(1.224744871391589*fr[18]+0.9486832980505137*fr[10]+0.5477225575051661*fr[5])+nuVtSq[1]*(0.9486832980505137*fr[17]+1.369306393762915*fr[14]+0.5477225575051661*fr[13]+1.060660171779821*fr[6]+0.6123724356957944*fr[3]); 
    incr2[11] = nuVtSq[1]*(1.224744871391589*fr[12]+0.9486832980505138*fr[4]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.060660171779821*fr[11]+0.6123724356957944*fr[7])+nuVtSq[2]*(0.6776309271789384*fr[11]+1.369306393762915*fr[8]+0.3912303982179757*fr[7]+1.060660171779821*fr[2]+0.6123724356957944*fr[0]); 
    incr2[12] = nuVtSq[0]*(5.303300858899105*fr[12]+4.107919181288746*fr[4]+2.371708245126284*fr[1])+nuVtSq[2]*(4.743416490252569*fr[12]+3.674234614174767*fr[4]+2.121320343559642*fr[1])+nuVtSq[1]*(3.674234614174766*fr[11]+5.303300858899106*fr[8]+2.121320343559642*fr[7]+4.107919181288746*fr[2]+2.371708245126284*fr[0]); 
    incr2[14] = nuVtSq[1]*(5.303300858899106*fr[18]+4.107919181288745*fr[10]+2.371708245126284*fr[5])+nuVtSq[2]*(4.107919181288745*fr[17]+2.371708245126284*fr[13])+nuVtSq[0]*(5.303300858899105*fr[14]+4.107919181288745*fr[6]+2.371708245126284*fr[3]); 
    incr2[16] = nuVtSq[1]*(1.060660171779821*fr[19]+0.6123724356957944*fr[15])+nuVtSq[0]*(1.060660171779821*fr[16]+0.6123724356957944*fr[9]); 
    incr2[17] = nuVtSq[1]*(1.224744871391589*fr[18]+0.9486832980505137*fr[10]+0.5477225575051661*fr[5])+nuVtSq[0]*(1.060660171779821*fr[17]+0.6123724356957944*fr[13])+nuVtSq[2]*(0.6776309271789384*fr[17]+1.369306393762915*fr[14]+0.3912303982179757*fr[13]+1.060660171779821*fr[6]+0.6123724356957944*fr[3]); 
    incr2[18] = nuVtSq[0]*(5.303300858899105*fr[18]+4.107919181288745*fr[10]+2.371708245126284*fr[5])+nuVtSq[2]*(4.743416490252569*fr[18]+3.674234614174766*fr[10]+2.121320343559642*fr[5])+nuVtSq[1]*(3.674234614174766*fr[17]+5.303300858899106*fr[14]+2.121320343559642*fr[13]+4.107919181288745*fr[6]+2.371708245126284*fr[3]); 
    incr2[19] = nuVtSq[0]*(1.060660171779821*fr[19]+0.6123724356957944*fr[15])+nuVtSq[2]*(0.9486832980505137*fr[19]+0.5477225575051661*fr[15])+nuVtSq[1]*(1.060660171779821*fr[16]+0.6123724356957944*fr[9]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[11] += incr2[11]*rdvSq4r; 
    outr[12] += incr2[12]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[17] += incr2[17]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 
    outr[19] += incr2[19]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[1]*((-1.369306393762915*fl[12])+1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[2]*(1.060660171779821*fl[11]-0.6123724356957944*fl[7])+nuVtSq[0]*((-1.369306393762915*fl[8])+1.060660171779821*fl[2]-0.6123724356957944*fl[0]); 
    incr2[4] = nuVtSq[2]*((-1.224744871391589*fl[12])+0.9486832980505137*fl[4]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[12])+1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[1]*(0.9486832980505138*fl[11]-1.369306393762915*fl[8]-0.5477225575051661*fl[7]+1.060660171779821*fl[2]-0.6123724356957944*fl[0]); 
    incr2[6] = nuVtSq[1]*((-1.369306393762915*fl[18])+1.060660171779821*fl[10]-0.6123724356957944*fl[5])+nuVtSq[2]*(1.060660171779821*fl[17]-0.6123724356957944*fl[13])+nuVtSq[0]*((-1.369306393762915*fl[14])+1.060660171779821*fl[6]-0.6123724356957944*fl[3]); 
    incr2[8] = nuVtSq[1]*(5.303300858899106*fl[12]-4.107919181288745*fl[4]+2.371708245126284*fl[1])+nuVtSq[2]*(2.371708245126284*fl[7]-4.107919181288745*fl[11])+nuVtSq[0]*(5.303300858899105*fl[8]-4.107919181288745*fl[2]+2.371708245126284*fl[0]); 
    incr2[10] = nuVtSq[2]*((-1.224744871391589*fl[18])+0.9486832980505137*fl[10]-0.5477225575051661*fl[5])+nuVtSq[0]*((-1.369306393762915*fl[18])+1.060660171779821*fl[10]-0.6123724356957944*fl[5])+nuVtSq[1]*(0.9486832980505137*fl[17]-1.369306393762915*fl[14]-0.5477225575051661*fl[13]+1.060660171779821*fl[6]-0.6123724356957944*fl[3]); 
    incr2[11] = nuVtSq[1]*((-1.224744871391589*fl[12])+0.9486832980505138*fl[4]-0.5477225575051661*fl[1])+nuVtSq[0]*(1.060660171779821*fl[11]-0.6123724356957944*fl[7])+nuVtSq[2]*(0.6776309271789384*fl[11]-1.369306393762915*fl[8]-0.3912303982179757*fl[7]+1.060660171779821*fl[2]-0.6123724356957944*fl[0]); 
    incr2[12] = nuVtSq[0]*(5.303300858899105*fl[12]-4.107919181288746*fl[4]+2.371708245126284*fl[1])+nuVtSq[2]*(4.743416490252569*fl[12]-3.674234614174767*fl[4]+2.121320343559642*fl[1])+nuVtSq[1]*((-3.674234614174766*fl[11])+5.303300858899106*fl[8]+2.121320343559642*fl[7]-4.107919181288746*fl[2]+2.371708245126284*fl[0]); 
    incr2[14] = nuVtSq[1]*(5.303300858899106*fl[18]-4.107919181288745*fl[10]+2.371708245126284*fl[5])+nuVtSq[2]*(2.371708245126284*fl[13]-4.107919181288745*fl[17])+nuVtSq[0]*(5.303300858899105*fl[14]-4.107919181288745*fl[6]+2.371708245126284*fl[3]); 
    incr2[16] = nuVtSq[1]*(1.060660171779821*fl[19]-0.6123724356957944*fl[15])+nuVtSq[0]*(1.060660171779821*fl[16]-0.6123724356957944*fl[9]); 
    incr2[17] = nuVtSq[1]*((-1.224744871391589*fl[18])+0.9486832980505137*fl[10]-0.5477225575051661*fl[5])+nuVtSq[0]*(1.060660171779821*fl[17]-0.6123724356957944*fl[13])+nuVtSq[2]*(0.6776309271789384*fl[17]-1.369306393762915*fl[14]-0.3912303982179757*fl[13]+1.060660171779821*fl[6]-0.6123724356957944*fl[3]); 
    incr2[18] = nuVtSq[0]*(5.303300858899105*fl[18]-4.107919181288745*fl[10]+2.371708245126284*fl[5])+nuVtSq[2]*(4.743416490252569*fl[18]-3.674234614174766*fl[10]+2.121320343559642*fl[5])+nuVtSq[1]*((-3.674234614174766*fl[17])+5.303300858899106*fl[14]+2.121320343559642*fl[13]-4.107919181288745*fl[6]+2.371708245126284*fl[3]); 
    incr2[19] = nuVtSq[0]*(1.060660171779821*fl[19]-0.6123724356957944*fl[15])+nuVtSq[2]*(0.9486832980505137*fl[19]-0.5477225575051661*fl[15])+nuVtSq[1]*(1.060660171779821*fl[16]-0.6123724356957944*fl[9]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[11] += incr2[11]*rdvSq4l; 
    outl[12] += incr2[12]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[17] += incr2[17]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 
    outl[19] += incr2[19]*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[8]; 
  if (idxr[2] == 1) {



  } else {



  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double mufac, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[20]; 
  if (idxr[2] == 1) {



  } else {



  }
  return 0.0; 
} 
