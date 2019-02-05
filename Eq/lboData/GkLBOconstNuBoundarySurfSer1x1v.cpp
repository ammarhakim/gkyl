#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu:        constant collisionality. 
  // vMuMidMax: maximum midpoint value of v-u. 
  // u[1*2]:    bulk velocity (in 1 directions). 
  // vtSq[2]:   thermal speed squared. 
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  if (idxr[1] == 1) {

    outr[2] += (vtSq[0]*(1.369306393762915*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3]))*rdvSq4nuR; 
    outr[3] += (vtSq[1]*(1.369306393762915*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3]))*rdvSq4nuR; 
    outr[4] += (vtSq[0]*((-5.303300858899105*fr[4])+4.107919181288745*fr[2]-2.371708245126284*fr[0])+vtSq[1]*(4.107919181288745*fr[3]-2.371708245126284*fr[1]))*rdvSq4nuR; 

  } else {

    outl[2] += (vtSq[0]*((-1.369306393762915*fl[4])-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[1]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1]))*rdvSq4nuL; 
    outl[3] += (vtSq[1]*((-1.369306393762915*fl[4])-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[0]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1]))*rdvSq4nuL; 
    outl[4] += (vtSq[0]*((-5.303300858899105*fl[4])-4.107919181288745*fl[2]-2.371708245126284*fl[0])+vtSq[1]*((-4.107919181288745*fl[3])-2.371708245126284*fl[1]))*rdvSq4nuL; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu:        constant collisionality. 
  // vMuMidMax: maximum midpoint value of v-u. 
  // u[1*3]:    bulk velocity (in 1 directions). 
  // vtSq[3]:   thermal speed squared. 
  // fl/fr:     Distribution function in left/right cells 
  // outl/outr: Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  if (idxr[1] == 1) {

    outr[2] += (vtSq[1]*(1.369306393762915*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[1])+vtSq[2]*(0.6123724356957944*fr[4]-1.060660171779821*fr[6])+vtSq[0]*(1.369306393762915*fr[5]-1.060660171779821*fr[2]+0.6123724356957944*fr[0]))*rdvSq4nuR; 
    outr[3] += (vtSq[0]*(1.369306393762915*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[1])+vtSq[2]*(1.224744871391589*fr[7]-0.9486832980505137*fr[3]+0.5477225575051661*fr[1])+vtSq[1]*((-0.9486832980505138*fr[6])+1.369306393762915*fr[5]+0.5477225575051661*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0]))*rdvSq4nuR; 
    outr[5] += (vtSq[1]*((-5.303300858899106*fr[7])+4.107919181288745*fr[3]-2.371708245126284*fr[1])+vtSq[2]*(4.107919181288745*fr[6]-2.371708245126284*fr[4])+vtSq[0]*((-5.303300858899105*fr[5])+4.107919181288745*fr[2]-2.371708245126284*fr[0]))*rdvSq4nuR; 
    outr[6] += (vtSq[1]*(1.224744871391589*fr[7]-0.9486832980505138*fr[3]+0.5477225575051661*fr[1])+vtSq[2]*((-0.6776309271789384*fr[6])+1.369306393762915*fr[5]+0.3912303982179757*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[0]*(0.6123724356957944*fr[4]-1.060660171779821*fr[6]))*rdvSq4nuR; 
    outr[7] += (vtSq[2]*((-4.743416490252569*fr[7])+3.674234614174767*fr[3]-2.121320343559642*fr[1])+vtSq[0]*((-5.303300858899105*fr[7])+4.107919181288746*fr[3]-2.371708245126284*fr[1])+vtSq[1]*(3.674234614174766*fr[6]-5.303300858899106*fr[5]-2.121320343559642*fr[4]+4.107919181288746*fr[2]-2.371708245126284*fr[0]))*rdvSq4nuR; 

  } else {

    outl[2] += (vtSq[1]*((-1.369306393762915*fl[7])-1.060660171779821*fl[3]-0.6123724356957944*fl[1])+vtSq[2]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[4])+vtSq[0]*((-1.369306393762915*fl[5])-1.060660171779821*fl[2]-0.6123724356957944*fl[0]))*rdvSq4nuL; 
    outl[3] += (vtSq[2]*((-1.224744871391589*fl[7])-0.9486832980505137*fl[3]-0.5477225575051661*fl[1])+vtSq[0]*((-1.369306393762915*fl[7])-1.060660171779821*fl[3]-0.6123724356957944*fl[1])+vtSq[1]*((-0.9486832980505138*fl[6])-1.369306393762915*fl[5]-0.5477225575051661*fl[4]-1.060660171779821*fl[2]-0.6123724356957944*fl[0]))*rdvSq4nuL; 
    outl[5] += (vtSq[1]*((-5.303300858899106*fl[7])-4.107919181288745*fl[3]-2.371708245126284*fl[1])+vtSq[2]*((-4.107919181288745*fl[6])-2.371708245126284*fl[4])+vtSq[0]*((-5.303300858899105*fl[5])-4.107919181288745*fl[2]-2.371708245126284*fl[0]))*rdvSq4nuL; 
    outl[6] += (vtSq[1]*((-1.224744871391589*fl[7])-0.9486832980505138*fl[3]-0.5477225575051661*fl[1])+vtSq[2]*((-0.6776309271789384*fl[6])-1.369306393762915*fl[5]-0.3912303982179757*fl[4]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[4]))*rdvSq4nuL; 
    outl[7] += (vtSq[2]*((-4.743416490252569*fl[7])-3.674234614174767*fl[3]-2.121320343559642*fl[1])+vtSq[0]*((-5.303300858899105*fl[7])-4.107919181288746*fl[3]-2.371708245126284*fl[1])+vtSq[1]*((-3.674234614174766*fl[6])-5.303300858899106*fl[5]-2.121320343559642*fl[4]-4.107919181288746*fl[2]-2.371708245126284*fl[0]))*rdvSq4nuL; 

  }
  return 0.0; 
} 
