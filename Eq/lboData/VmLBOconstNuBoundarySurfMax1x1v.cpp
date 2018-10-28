#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf1x1vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates.
  // dxv[2]:     Cell spacing.
  // idx[2]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[2]:       bulk velocity (in 1 directions). 
  // vtSq[2]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double incr2[3]; 
  if (idxr[1] == 1) {

    incr2[2] = vtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2])+0.6123724356957944*fr[1]*vtSq[1]; 

    outr[2] += incr2[2]*rdvSq4nuR; 

  } else {

    incr2[2] = vtSq[0]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*vtSq[1]; 

    outl[2] += incr2[2]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x1vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates.
  // dxv[2]:     Cell spacing.
  // idx[2]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[3]:       bulk velocity (in 1 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double incr2[6]; 
  if (idxr[1] == 1) {

    incr2[2] = vtSq[0]*(1.369306393762915*fr[5]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*vtSq[2]*fr[4]+vtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3]); 
    incr2[3] = vtSq[1]*(1.369306393762915*fr[5]+0.5477225575051661*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[3])+vtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3]); 
    incr2[5] = vtSq[0]*((-5.303300858899105*fr[5])+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*vtSq[2]*fr[4]+vtSq[1]*(4.107919181288745*fr[3]-2.371708245126284*fr[1]); 

    outr[2] += incr2[2]*rdvSq4nuR; 
    outr[3] += incr2[3]*rdvSq4nuR; 
    outr[5] += incr2[5]*rdvSq4nuR; 

  } else {

    incr2[2] = vtSq[0]*((-1.369306393762915*fl[5])-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*vtSq[2]*fl[4]+vtSq[1]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1]); 
    incr2[3] = vtSq[1]*((-1.369306393762915*fl[5])-0.5477225575051661*fl[4]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[2]*((-0.9486832980505137*fl[3])-0.5477225575051661*fl[1])+vtSq[0]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1]); 
    incr2[5] = vtSq[0]*((-5.303300858899105*fl[5])-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*vtSq[2]*fl[4]+vtSq[1]*((-4.107919181288745*fl[3])-2.371708245126284*fl[1]); 

    outl[2] += incr2[2]*rdvSq4nuL; 
    outl[3] += incr2[3]*rdvSq4nuL; 
    outl[5] += incr2[5]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x1vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates.
  // dxv[2]:     Cell spacing.
  // idx[2]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[4]:       bulk velocity (in 1 directions). 
  // vtSq[4]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[1]*dxvl[1]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[1]*dxvr[1]); 

  double incr2[10]; 
  if (idxr[1] == 1) {

    incr2[2] = vtSq[0]*((-1.620185174601965*fr[9])+1.369306393762915*fr[5]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*vtSq[3]*fr[8]+vtSq[1]*(1.369306393762915*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[1])+vtSq[2]*(0.6123724356957944*fr[4]-1.060660171779821*fr[6]); 
    incr2[3] = vtSq[1]*((-1.620185174601965*fr[9])-0.9486832980505138*fr[6]+1.369306393762915*fr[5]+0.5477225575051661*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[2]*(0.537852874200477*fr[8]+1.224744871391589*fr[7]-0.9486832980505137*fr[3]+0.5477225575051661*fr[1])+vtSq[0]*(1.369306393762915*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[1])+vtSq[3]*(0.537852874200477*fr[4]-0.931588505112178*fr[6]); 
    incr2[5] = vtSq[0]*(6.274950199005565*fr[9]-5.303300858899105*fr[5]+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*vtSq[3]*fr[8]+vtSq[1]*((-5.303300858899106*fr[7])+4.107919181288745*fr[3]-2.371708245126284*fr[1])+vtSq[2]*(4.107919181288745*fr[6]-2.371708245126284*fr[4]); 
    incr2[6] = vtSq[2]*((-1.620185174601965*fr[9])-0.6776309271789384*fr[6]+1.369306393762915*fr[5]+0.3912303982179757*fr[4]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+vtSq[1]*(0.5378528742004769*fr[8]+1.224744871391589*fr[7]-0.9486832980505138*fr[3]+0.5477225575051661*fr[1])+vtSq[3]*(0.3651483716701107*fr[8]+1.202675588605909*fr[7]-0.931588505112178*fr[3]+0.5378528742004769*fr[1])+vtSq[0]*(0.6123724356957944*fr[4]-1.060660171779821*fr[6]); 
    incr2[7] = vtSq[1]*(6.274950199005565*fr[9]+3.674234614174766*fr[6]-5.303300858899106*fr[5]-2.121320343559642*fr[4]+4.107919181288746*fr[2]-2.371708245126284*fr[0])+vtSq[2]*((-2.08309522448824*fr[8])-4.743416490252569*fr[7]+3.674234614174767*fr[3]-2.121320343559642*fr[1])+vtSq[0]*((-5.303300858899105*fr[7])+4.107919181288746*fr[3]-2.371708245126284*fr[1])+vtSq[3]*(3.608026765817728*fr[6]-2.08309522448824*fr[4]); 
    incr2[9] = vtSq[0]*((-14.8492424049175*fr[9])+12.54990039801113*fr[5]-9.721111047611789*fr[2]+5.612486080160912*fr[0])+5.612486080160912*vtSq[3]*fr[8]+vtSq[1]*(12.54990039801113*fr[7]-9.721111047611789*fr[3]+5.612486080160912*fr[1])+vtSq[2]*(5.612486080160912*fr[4]-9.721111047611789*fr[6]); 

    outr[2] += incr2[2]*rdvSq4nuR; 
    outr[3] += incr2[3]*rdvSq4nuR; 
    outr[5] += incr2[5]*rdvSq4nuR; 
    outr[6] += incr2[6]*rdvSq4nuR; 
    outr[7] += incr2[7]*rdvSq4nuR; 
    outr[9] += incr2[9]*rdvSq4nuR; 

  } else {

    incr2[2] = vtSq[0]*((-1.620185174601965*fl[9])-1.369306393762915*fl[5]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*vtSq[3]*fl[8]+vtSq[1]*((-1.369306393762915*fl[7])-1.060660171779821*fl[3]-0.6123724356957944*fl[1])+vtSq[2]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[4]); 
    incr2[3] = vtSq[1]*((-1.620185174601965*fl[9])-0.9486832980505138*fl[6]-1.369306393762915*fl[5]-0.5477225575051661*fl[4]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[2]*((-0.537852874200477*fl[8])-1.224744871391589*fl[7]-0.9486832980505137*fl[3]-0.5477225575051661*fl[1])+vtSq[0]*((-1.369306393762915*fl[7])-1.060660171779821*fl[3]-0.6123724356957944*fl[1])+vtSq[3]*((-0.931588505112178*fl[6])-0.537852874200477*fl[4]); 
    incr2[5] = vtSq[0]*((-6.274950199005565*fl[9])-5.303300858899105*fl[5]-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*vtSq[3]*fl[8]+vtSq[1]*((-5.303300858899106*fl[7])-4.107919181288745*fl[3]-2.371708245126284*fl[1])+vtSq[2]*((-4.107919181288745*fl[6])-2.371708245126284*fl[4]); 
    incr2[6] = vtSq[2]*((-1.620185174601965*fl[9])-0.6776309271789384*fl[6]-1.369306393762915*fl[5]-0.3912303982179757*fl[4]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+vtSq[3]*((-0.3651483716701107*fl[8])-1.202675588605909*fl[7]-0.931588505112178*fl[3]-0.5378528742004769*fl[1])+vtSq[1]*((-0.5378528742004769*fl[8])-1.224744871391589*fl[7]-0.9486832980505138*fl[3]-0.5477225575051661*fl[1])+vtSq[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[4]); 
    incr2[7] = vtSq[1]*((-6.274950199005565*fl[9])-3.674234614174766*fl[6]-5.303300858899106*fl[5]-2.121320343559642*fl[4]-4.107919181288746*fl[2]-2.371708245126284*fl[0])+vtSq[2]*((-2.08309522448824*fl[8])-4.743416490252569*fl[7]-3.674234614174767*fl[3]-2.121320343559642*fl[1])+vtSq[0]*((-5.303300858899105*fl[7])-4.107919181288746*fl[3]-2.371708245126284*fl[1])+vtSq[3]*((-3.608026765817728*fl[6])-2.08309522448824*fl[4]); 
    incr2[9] = vtSq[0]*((-14.8492424049175*fl[9])-12.54990039801113*fl[5]-9.721111047611789*fl[2]-5.612486080160912*fl[0])-5.612486080160912*vtSq[3]*fl[8]+vtSq[1]*((-12.54990039801113*fl[7])-9.721111047611789*fl[3]-5.612486080160912*fl[1])+vtSq[2]*((-9.721111047611789*fl[6])-5.612486080160912*fl[4]); 

    outl[2] += incr2[2]*rdvSq4nuL; 
    outl[3] += incr2[3]*rdvSq4nuL; 
    outl[5] += incr2[5]*rdvSq4nuL; 
    outl[6] += incr2[6]*rdvSq4nuL; 
    outl[7] += incr2[7]*rdvSq4nuL; 
    outl[9] += incr2[9]*rdvSq4nuL; 

  }
  return 0.0; 
} 
