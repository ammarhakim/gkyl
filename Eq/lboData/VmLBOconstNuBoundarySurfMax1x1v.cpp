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

  if (idxr[1] == 1) {

    outr[2] += (-1.060660171779821*vtSq[0]*fr[2]*rdvSq4nuR)+0.6123724356957944*fr[1]*vtSq[1]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[2] += (-1.060660171779821*vtSq[0]*fl[2]*rdvSq4nuL)-0.6123724356957944*fl[1]*vtSq[1]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[0]*rdvSq4nuL; 

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

  if (idxr[1] == 1) {

    outr[2] += 1.369306393762915*vtSq[0]*fr[5]*rdvSq4nuR+0.6123724356957944*vtSq[2]*fr[4]*rdvSq4nuR-1.060660171779821*vtSq[1]*fr[3]*rdvSq4nuR-1.060660171779821*vtSq[0]*fr[2]*rdvSq4nuR+0.6123724356957944*fr[1]*vtSq[1]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[3] += 1.369306393762915*vtSq[1]*fr[5]*rdvSq4nuR+0.5477225575051661*vtSq[1]*fr[4]*rdvSq4nuR-0.9486832980505137*vtSq[2]*fr[3]*rdvSq4nuR-1.060660171779821*vtSq[0]*fr[3]*rdvSq4nuR+0.5477225575051661*fr[1]*vtSq[2]*rdvSq4nuR-1.060660171779821*vtSq[1]*fr[2]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[1]*rdvSq4nuR+0.6123724356957944*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[5] += (-5.303300858899105*vtSq[0]*fr[5]*rdvSq4nuR)-2.371708245126284*vtSq[2]*fr[4]*rdvSq4nuR+4.107919181288745*vtSq[1]*fr[3]*rdvSq4nuR+4.107919181288745*vtSq[0]*fr[2]*rdvSq4nuR-2.371708245126284*fr[1]*vtSq[1]*rdvSq4nuR-2.371708245126284*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[2] += (-1.369306393762915*vtSq[0]*fl[5]*rdvSq4nuL)-0.6123724356957944*vtSq[2]*fl[4]*rdvSq4nuL-1.060660171779821*vtSq[1]*fl[3]*rdvSq4nuL-1.060660171779821*vtSq[0]*fl[2]*rdvSq4nuL-0.6123724356957944*fl[1]*vtSq[1]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[3] += (-1.369306393762915*vtSq[1]*fl[5]*rdvSq4nuL)-0.5477225575051661*vtSq[1]*fl[4]*rdvSq4nuL-0.9486832980505137*vtSq[2]*fl[3]*rdvSq4nuL-1.060660171779821*vtSq[0]*fl[3]*rdvSq4nuL-0.5477225575051661*fl[1]*vtSq[2]*rdvSq4nuL-1.060660171779821*vtSq[1]*fl[2]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[1]*rdvSq4nuL-0.6123724356957944*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[5] += (-5.303300858899105*vtSq[0]*fl[5]*rdvSq4nuL)-2.371708245126284*vtSq[2]*fl[4]*rdvSq4nuL-4.107919181288745*vtSq[1]*fl[3]*rdvSq4nuL-4.107919181288745*vtSq[0]*fl[2]*rdvSq4nuL-2.371708245126284*fl[1]*vtSq[1]*rdvSq4nuL-2.371708245126284*fl[0]*vtSq[0]*rdvSq4nuL; 

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

  if (idxr[1] == 1) {

    outr[2] += (-1.620185174601965*vtSq[0]*fr[9]*rdvSq4nuR)+0.6123724356957944*vtSq[3]*fr[8]*rdvSq4nuR+1.369306393762915*vtSq[1]*fr[7]*rdvSq4nuR-1.060660171779821*vtSq[2]*fr[6]*rdvSq4nuR+1.369306393762915*vtSq[0]*fr[5]*rdvSq4nuR+0.6123724356957944*vtSq[2]*fr[4]*rdvSq4nuR-1.060660171779821*vtSq[1]*fr[3]*rdvSq4nuR-1.060660171779821*vtSq[0]*fr[2]*rdvSq4nuR+0.6123724356957944*fr[1]*vtSq[1]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[3] += (-1.620185174601965*vtSq[1]*fr[9]*rdvSq4nuR)+0.537852874200477*vtSq[2]*fr[8]*rdvSq4nuR+1.224744871391589*vtSq[2]*fr[7]*rdvSq4nuR+1.369306393762915*vtSq[0]*fr[7]*rdvSq4nuR-0.931588505112178*vtSq[3]*fr[6]*rdvSq4nuR-0.9486832980505138*vtSq[1]*fr[6]*rdvSq4nuR+1.369306393762915*vtSq[1]*fr[5]*rdvSq4nuR+0.537852874200477*vtSq[3]*fr[4]*rdvSq4nuR+0.5477225575051661*vtSq[1]*fr[4]*rdvSq4nuR-0.9486832980505137*vtSq[2]*fr[3]*rdvSq4nuR-1.060660171779821*vtSq[0]*fr[3]*rdvSq4nuR+0.5477225575051661*fr[1]*vtSq[2]*rdvSq4nuR-1.060660171779821*vtSq[1]*fr[2]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[1]*rdvSq4nuR+0.6123724356957944*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[5] += 6.274950199005565*vtSq[0]*fr[9]*rdvSq4nuR-2.371708245126284*vtSq[3]*fr[8]*rdvSq4nuR-5.303300858899106*vtSq[1]*fr[7]*rdvSq4nuR+4.107919181288745*vtSq[2]*fr[6]*rdvSq4nuR-5.303300858899105*vtSq[0]*fr[5]*rdvSq4nuR-2.371708245126284*vtSq[2]*fr[4]*rdvSq4nuR+4.107919181288745*vtSq[1]*fr[3]*rdvSq4nuR+4.107919181288745*vtSq[0]*fr[2]*rdvSq4nuR-2.371708245126284*fr[1]*vtSq[1]*rdvSq4nuR-2.371708245126284*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[6] += (-1.620185174601965*vtSq[2]*fr[9]*rdvSq4nuR)+0.3651483716701107*vtSq[3]*fr[8]*rdvSq4nuR+0.5378528742004769*vtSq[1]*fr[8]*rdvSq4nuR+1.202675588605909*vtSq[3]*fr[7]*rdvSq4nuR+1.224744871391589*vtSq[1]*fr[7]*rdvSq4nuR-0.6776309271789384*vtSq[2]*fr[6]*rdvSq4nuR-1.060660171779821*vtSq[0]*fr[6]*rdvSq4nuR+1.369306393762915*vtSq[2]*fr[5]*rdvSq4nuR+0.3912303982179757*vtSq[2]*fr[4]*rdvSq4nuR+0.6123724356957944*vtSq[0]*fr[4]*rdvSq4nuR-0.931588505112178*fr[3]*vtSq[3]*rdvSq4nuR+0.5378528742004769*fr[1]*vtSq[3]*rdvSq4nuR-0.9486832980505138*vtSq[1]*fr[3]*rdvSq4nuR-1.060660171779821*fr[2]*vtSq[2]*rdvSq4nuR+0.6123724356957944*fr[0]*vtSq[2]*rdvSq4nuR+0.5477225575051661*fr[1]*vtSq[1]*rdvSq4nuR; 
    outr[7] += 6.274950199005565*vtSq[1]*fr[9]*rdvSq4nuR-2.08309522448824*vtSq[2]*fr[8]*rdvSq4nuR-4.743416490252569*vtSq[2]*fr[7]*rdvSq4nuR-5.303300858899105*vtSq[0]*fr[7]*rdvSq4nuR+3.608026765817728*vtSq[3]*fr[6]*rdvSq4nuR+3.674234614174766*vtSq[1]*fr[6]*rdvSq4nuR-5.303300858899106*vtSq[1]*fr[5]*rdvSq4nuR-2.08309522448824*vtSq[3]*fr[4]*rdvSq4nuR-2.121320343559642*vtSq[1]*fr[4]*rdvSq4nuR+3.674234614174767*vtSq[2]*fr[3]*rdvSq4nuR+4.107919181288746*vtSq[0]*fr[3]*rdvSq4nuR-2.121320343559642*fr[1]*vtSq[2]*rdvSq4nuR+4.107919181288746*vtSq[1]*fr[2]*rdvSq4nuR-2.371708245126284*fr[0]*vtSq[1]*rdvSq4nuR-2.371708245126284*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[9] += (-14.8492424049175*vtSq[0]*fr[9]*rdvSq4nuR)+5.612486080160912*vtSq[3]*fr[8]*rdvSq4nuR+12.54990039801113*vtSq[1]*fr[7]*rdvSq4nuR-9.721111047611789*vtSq[2]*fr[6]*rdvSq4nuR+12.54990039801113*vtSq[0]*fr[5]*rdvSq4nuR+5.612486080160912*vtSq[2]*fr[4]*rdvSq4nuR-9.721111047611789*vtSq[1]*fr[3]*rdvSq4nuR-9.721111047611789*vtSq[0]*fr[2]*rdvSq4nuR+5.612486080160912*fr[1]*vtSq[1]*rdvSq4nuR+5.612486080160912*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[2] += (-1.620185174601965*vtSq[0]*fl[9]*rdvSq4nuL)-0.6123724356957944*vtSq[3]*fl[8]*rdvSq4nuL-1.369306393762915*vtSq[1]*fl[7]*rdvSq4nuL-1.060660171779821*vtSq[2]*fl[6]*rdvSq4nuL-1.369306393762915*vtSq[0]*fl[5]*rdvSq4nuL-0.6123724356957944*vtSq[2]*fl[4]*rdvSq4nuL-1.060660171779821*vtSq[1]*fl[3]*rdvSq4nuL-1.060660171779821*vtSq[0]*fl[2]*rdvSq4nuL-0.6123724356957944*fl[1]*vtSq[1]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[3] += (-1.620185174601965*vtSq[1]*fl[9]*rdvSq4nuL)-0.537852874200477*vtSq[2]*fl[8]*rdvSq4nuL-1.224744871391589*vtSq[2]*fl[7]*rdvSq4nuL-1.369306393762915*vtSq[0]*fl[7]*rdvSq4nuL-0.931588505112178*vtSq[3]*fl[6]*rdvSq4nuL-0.9486832980505138*vtSq[1]*fl[6]*rdvSq4nuL-1.369306393762915*vtSq[1]*fl[5]*rdvSq4nuL-0.537852874200477*vtSq[3]*fl[4]*rdvSq4nuL-0.5477225575051661*vtSq[1]*fl[4]*rdvSq4nuL-0.9486832980505137*vtSq[2]*fl[3]*rdvSq4nuL-1.060660171779821*vtSq[0]*fl[3]*rdvSq4nuL-0.5477225575051661*fl[1]*vtSq[2]*rdvSq4nuL-1.060660171779821*vtSq[1]*fl[2]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[1]*rdvSq4nuL-0.6123724356957944*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[5] += (-6.274950199005565*vtSq[0]*fl[9]*rdvSq4nuL)-2.371708245126284*vtSq[3]*fl[8]*rdvSq4nuL-5.303300858899106*vtSq[1]*fl[7]*rdvSq4nuL-4.107919181288745*vtSq[2]*fl[6]*rdvSq4nuL-5.303300858899105*vtSq[0]*fl[5]*rdvSq4nuL-2.371708245126284*vtSq[2]*fl[4]*rdvSq4nuL-4.107919181288745*vtSq[1]*fl[3]*rdvSq4nuL-4.107919181288745*vtSq[0]*fl[2]*rdvSq4nuL-2.371708245126284*fl[1]*vtSq[1]*rdvSq4nuL-2.371708245126284*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[6] += (-1.620185174601965*vtSq[2]*fl[9]*rdvSq4nuL)-0.3651483716701107*vtSq[3]*fl[8]*rdvSq4nuL-0.5378528742004769*vtSq[1]*fl[8]*rdvSq4nuL-1.202675588605909*vtSq[3]*fl[7]*rdvSq4nuL-1.224744871391589*vtSq[1]*fl[7]*rdvSq4nuL-0.6776309271789384*vtSq[2]*fl[6]*rdvSq4nuL-1.060660171779821*vtSq[0]*fl[6]*rdvSq4nuL-1.369306393762915*vtSq[2]*fl[5]*rdvSq4nuL-0.3912303982179757*vtSq[2]*fl[4]*rdvSq4nuL-0.6123724356957944*vtSq[0]*fl[4]*rdvSq4nuL-0.931588505112178*fl[3]*vtSq[3]*rdvSq4nuL-0.5378528742004769*fl[1]*vtSq[3]*rdvSq4nuL-0.9486832980505138*vtSq[1]*fl[3]*rdvSq4nuL-1.060660171779821*fl[2]*vtSq[2]*rdvSq4nuL-0.6123724356957944*fl[0]*vtSq[2]*rdvSq4nuL-0.5477225575051661*fl[1]*vtSq[1]*rdvSq4nuL; 
    outl[7] += (-6.274950199005565*vtSq[1]*fl[9]*rdvSq4nuL)-2.08309522448824*vtSq[2]*fl[8]*rdvSq4nuL-4.743416490252569*vtSq[2]*fl[7]*rdvSq4nuL-5.303300858899105*vtSq[0]*fl[7]*rdvSq4nuL-3.608026765817728*vtSq[3]*fl[6]*rdvSq4nuL-3.674234614174766*vtSq[1]*fl[6]*rdvSq4nuL-5.303300858899106*vtSq[1]*fl[5]*rdvSq4nuL-2.08309522448824*vtSq[3]*fl[4]*rdvSq4nuL-2.121320343559642*vtSq[1]*fl[4]*rdvSq4nuL-3.674234614174767*vtSq[2]*fl[3]*rdvSq4nuL-4.107919181288746*vtSq[0]*fl[3]*rdvSq4nuL-2.121320343559642*fl[1]*vtSq[2]*rdvSq4nuL-4.107919181288746*vtSq[1]*fl[2]*rdvSq4nuL-2.371708245126284*fl[0]*vtSq[1]*rdvSq4nuL-2.371708245126284*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[9] += (-14.8492424049175*vtSq[0]*fl[9]*rdvSq4nuL)-5.612486080160912*vtSq[3]*fl[8]*rdvSq4nuL-12.54990039801113*vtSq[1]*fl[7]*rdvSq4nuL-9.721111047611789*vtSq[2]*fl[6]*rdvSq4nuL-12.54990039801113*vtSq[0]*fl[5]*rdvSq4nuL-5.612486080160912*vtSq[2]*fl[4]*rdvSq4nuL-9.721111047611789*vtSq[1]*fl[3]*rdvSq4nuL-9.721111047611789*vtSq[0]*fl[2]*rdvSq4nuL-5.612486080160912*fl[1]*vtSq[1]*rdvSq4nuL-5.612486080160912*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
