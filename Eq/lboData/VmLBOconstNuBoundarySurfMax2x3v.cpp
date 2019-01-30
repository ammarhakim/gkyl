#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf2x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += (-0.75*vtSq[0]*fr[3]*rdvSq4nuR)+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[3] += (-0.75*vtSq[0]*fl[3]*rdvSq4nuL)-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[2]*dxvl[2]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += 0.9682458365518543*vtSq[0]*fr[18]*rdvSq4nuR+0.4330127018922193*vtSq[5]*fr[17]*rdvSq4nuR+0.4330127018922193*vtSq[4]*fr[16]*rdvSq4nuR-0.75*vtSq[2]*fr[8]*rdvSq4nuR-0.75*vtSq[1]*fr[7]*rdvSq4nuR+0.4330127018922193*vtSq[3]*fr[6]*rdvSq4nuR-0.75*vtSq[0]*fr[3]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[7] += 0.9682458365518543*vtSq[1]*fr[18]*rdvSq4nuR+0.3872983346207416*vtSq[1]*fr[16]*rdvSq4nuR-0.75*vtSq[3]*fr[8]*rdvSq4nuR-0.6708203932499369*vtSq[4]*fr[7]*rdvSq4nuR-0.75*vtSq[0]*fr[7]*rdvSq4nuR+0.4330127018922193*vtSq[2]*fr[6]*rdvSq4nuR+0.3872983346207416*fr[1]*vtSq[4]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[3]*rdvSq4nuR-0.75*vtSq[1]*fr[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[1]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[8] += 0.9682458365518543*vtSq[2]*fr[18]*rdvSq4nuR+0.3872983346207416*vtSq[2]*fr[17]*rdvSq4nuR-0.6708203932499369*vtSq[5]*fr[8]*rdvSq4nuR-0.75*vtSq[0]*fr[8]*rdvSq4nuR-0.75*vtSq[3]*fr[7]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[6]*rdvSq4nuR+0.3872983346207416*fr[2]*vtSq[5]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[3]*rdvSq4nuR-0.75*vtSq[2]*fr[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[2]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[2]*rdvSq4nuR; 
    outr[11] += (-0.75*vtSq[0]*fr[11]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[10]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[9]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[4]*rdvSq4nuR; 
    outr[14] += (-0.75*vtSq[0]*fr[14]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[13]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[12]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[5]*rdvSq4nuR; 
    outr[18] += (-3.75*vtSq[0]*fr[18]*rdvSq4nuR)-1.677050983124842*vtSq[5]*fr[17]*rdvSq4nuR-1.677050983124842*vtSq[4]*fr[16]*rdvSq4nuR+2.904737509655563*vtSq[2]*fr[8]*rdvSq4nuR+2.904737509655563*vtSq[1]*fr[7]*rdvSq4nuR-1.677050983124842*vtSq[3]*fr[6]*rdvSq4nuR+2.904737509655563*vtSq[0]*fr[3]*rdvSq4nuR-1.677050983124842*fr[2]*vtSq[2]*rdvSq4nuR-1.677050983124842*fr[1]*vtSq[1]*rdvSq4nuR-1.677050983124842*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[3] += (-0.9682458365518543*vtSq[0]*fl[18]*rdvSq4nuL)-0.4330127018922193*vtSq[5]*fl[17]*rdvSq4nuL-0.4330127018922193*vtSq[4]*fl[16]*rdvSq4nuL-0.75*vtSq[2]*fl[8]*rdvSq4nuL-0.75*vtSq[1]*fl[7]*rdvSq4nuL-0.4330127018922193*vtSq[3]*fl[6]*rdvSq4nuL-0.75*vtSq[0]*fl[3]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[7] += (-0.9682458365518543*vtSq[1]*fl[18]*rdvSq4nuL)-0.3872983346207416*vtSq[1]*fl[16]*rdvSq4nuL-0.75*vtSq[3]*fl[8]*rdvSq4nuL-0.6708203932499369*vtSq[4]*fl[7]*rdvSq4nuL-0.75*vtSq[0]*fl[7]*rdvSq4nuL-0.4330127018922193*vtSq[2]*fl[6]*rdvSq4nuL-0.3872983346207416*fl[1]*vtSq[4]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[3]*rdvSq4nuL-0.75*vtSq[1]*fl[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[1]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[8] += (-0.9682458365518543*vtSq[2]*fl[18]*rdvSq4nuL)-0.3872983346207416*vtSq[2]*fl[17]*rdvSq4nuL-0.6708203932499369*vtSq[5]*fl[8]*rdvSq4nuL-0.75*vtSq[0]*fl[8]*rdvSq4nuL-0.75*vtSq[3]*fl[7]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[6]*rdvSq4nuL-0.3872983346207416*fl[2]*vtSq[5]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[3]*rdvSq4nuL-0.75*vtSq[2]*fl[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[2]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[2]*rdvSq4nuL; 
    outl[11] += (-0.75*vtSq[0]*fl[11]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[10]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[9]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[4]*rdvSq4nuL; 
    outl[14] += (-0.75*vtSq[0]*fl[14]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[13]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[12]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[5]*rdvSq4nuL; 
    outl[18] += (-3.75*vtSq[0]*fl[18]*rdvSq4nuL)-1.677050983124842*vtSq[5]*fl[17]*rdvSq4nuL-1.677050983124842*vtSq[4]*fl[16]*rdvSq4nuL-2.904737509655563*vtSq[2]*fl[8]*rdvSq4nuL-2.904737509655563*vtSq[1]*fl[7]*rdvSq4nuL-1.677050983124842*vtSq[3]*fl[6]*rdvSq4nuL-2.904737509655563*vtSq[0]*fl[3]*rdvSq4nuL-1.677050983124842*fl[2]*vtSq[2]*rdvSq4nuL-1.677050983124842*fl[1]*vtSq[1]*rdvSq4nuL-1.677050983124842*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[3]*dxvl[3]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += (-0.75*vtSq[0]*fr[4]*rdvSq4nuR)+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[4] += (-0.75*vtSq[0]*fl[4]*rdvSq4nuL)-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[3]*dxvl[3]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += 0.9682458365518543*vtSq[0]*fr[19]*rdvSq4nuR+0.4330127018922193*vtSq[5]*fr[17]*rdvSq4nuR+0.4330127018922193*vtSq[4]*fr[16]*rdvSq4nuR-0.75*vtSq[2]*fr[10]*rdvSq4nuR-0.75*vtSq[1]*fr[9]*rdvSq4nuR+0.4330127018922193*vtSq[3]*fr[6]*rdvSq4nuR-0.75*vtSq[0]*fr[4]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[9] += 0.9682458365518543*vtSq[1]*fr[19]*rdvSq4nuR+0.3872983346207416*vtSq[1]*fr[16]*rdvSq4nuR-0.75*vtSq[3]*fr[10]*rdvSq4nuR-0.6708203932499369*vtSq[4]*fr[9]*rdvSq4nuR-0.75*vtSq[0]*fr[9]*rdvSq4nuR+0.4330127018922193*vtSq[2]*fr[6]*rdvSq4nuR+0.3872983346207416*fr[1]*vtSq[4]*rdvSq4nuR-0.75*vtSq[1]*fr[4]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[1]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[10] += 0.9682458365518543*vtSq[2]*fr[19]*rdvSq4nuR+0.3872983346207416*vtSq[2]*fr[17]*rdvSq4nuR-0.6708203932499369*vtSq[5]*fr[10]*rdvSq4nuR-0.75*vtSq[0]*fr[10]*rdvSq4nuR-0.75*vtSq[3]*fr[9]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[6]*rdvSq4nuR+0.3872983346207416*fr[2]*vtSq[5]*rdvSq4nuR-0.75*vtSq[2]*fr[4]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[2]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[2]*rdvSq4nuR; 
    outr[11] += (-0.75*vtSq[0]*fr[11]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[8]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[7]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[3]*rdvSq4nuR; 
    outr[15] += (-0.75*vtSq[0]*fr[15]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[13]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[12]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[5]*rdvSq4nuR; 
    outr[19] += (-3.75*vtSq[0]*fr[19]*rdvSq4nuR)-1.677050983124842*vtSq[5]*fr[17]*rdvSq4nuR-1.677050983124842*vtSq[4]*fr[16]*rdvSq4nuR+2.904737509655563*vtSq[2]*fr[10]*rdvSq4nuR+2.904737509655563*vtSq[1]*fr[9]*rdvSq4nuR-1.677050983124842*vtSq[3]*fr[6]*rdvSq4nuR+2.904737509655563*vtSq[0]*fr[4]*rdvSq4nuR-1.677050983124842*fr[2]*vtSq[2]*rdvSq4nuR-1.677050983124842*fr[1]*vtSq[1]*rdvSq4nuR-1.677050983124842*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[4] += (-0.9682458365518543*vtSq[0]*fl[19]*rdvSq4nuL)-0.4330127018922193*vtSq[5]*fl[17]*rdvSq4nuL-0.4330127018922193*vtSq[4]*fl[16]*rdvSq4nuL-0.75*vtSq[2]*fl[10]*rdvSq4nuL-0.75*vtSq[1]*fl[9]*rdvSq4nuL-0.4330127018922193*vtSq[3]*fl[6]*rdvSq4nuL-0.75*vtSq[0]*fl[4]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[9] += (-0.9682458365518543*vtSq[1]*fl[19]*rdvSq4nuL)-0.3872983346207416*vtSq[1]*fl[16]*rdvSq4nuL-0.75*vtSq[3]*fl[10]*rdvSq4nuL-0.6708203932499369*vtSq[4]*fl[9]*rdvSq4nuL-0.75*vtSq[0]*fl[9]*rdvSq4nuL-0.4330127018922193*vtSq[2]*fl[6]*rdvSq4nuL-0.3872983346207416*fl[1]*vtSq[4]*rdvSq4nuL-0.75*vtSq[1]*fl[4]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[1]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[10] += (-0.9682458365518543*vtSq[2]*fl[19]*rdvSq4nuL)-0.3872983346207416*vtSq[2]*fl[17]*rdvSq4nuL-0.6708203932499369*vtSq[5]*fl[10]*rdvSq4nuL-0.75*vtSq[0]*fl[10]*rdvSq4nuL-0.75*vtSq[3]*fl[9]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[6]*rdvSq4nuL-0.3872983346207416*fl[2]*vtSq[5]*rdvSq4nuL-0.75*vtSq[2]*fl[4]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[2]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[2]*rdvSq4nuL; 
    outl[11] += (-0.75*vtSq[0]*fl[11]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[8]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[7]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[3]*rdvSq4nuL; 
    outl[15] += (-0.75*vtSq[0]*fl[15]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[13]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[12]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[5]*rdvSq4nuL; 
    outl[19] += (-3.75*vtSq[0]*fl[19]*rdvSq4nuL)-1.677050983124842*vtSq[5]*fl[17]*rdvSq4nuL-1.677050983124842*vtSq[4]*fl[16]*rdvSq4nuL-2.904737509655563*vtSq[2]*fl[10]*rdvSq4nuL-2.904737509655563*vtSq[1]*fl[9]*rdvSq4nuL-1.677050983124842*vtSq[3]*fl[6]*rdvSq4nuL-2.904737509655563*vtSq[0]*fl[4]*rdvSq4nuL-1.677050983124842*fl[2]*vtSq[2]*rdvSq4nuL-1.677050983124842*fl[1]*vtSq[1]*rdvSq4nuL-1.677050983124842*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[9]:       bulk velocity (in 3 directions). 
  // vtSq[3]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[4]*dxvl[4]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[4]*dxvr[4]); 

  if (idxr[4] == 1) {

    outr[5] += (-0.75*vtSq[0]*fr[5]*rdvSq4nuR)+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[5] += (-0.75*vtSq[0]*fl[5]*rdvSq4nuL)-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *u, const double *vtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:       Cell-center coordinates.
  // dxv[5]:     Cell spacing.
  // idx[5]:     current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // u[18]:       bulk velocity (in 3 directions). 
  // vtSq[6]:    thermal speed squared. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4nuL = 4.0*nu/(dxvl[4]*dxvl[4]); 
  double rdvSq4nuR = 4.0*nu/(dxvr[4]*dxvr[4]); 

  if (idxr[4] == 1) {

    outr[5] += 0.9682458365518543*vtSq[0]*fr[20]*rdvSq4nuR+0.4330127018922193*vtSq[5]*fr[17]*rdvSq4nuR+0.4330127018922193*vtSq[4]*fr[16]*rdvSq4nuR-0.75*vtSq[2]*fr[13]*rdvSq4nuR-0.75*vtSq[1]*fr[12]*rdvSq4nuR+0.4330127018922193*vtSq[3]*fr[6]*rdvSq4nuR-0.75*vtSq[0]*fr[5]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[2]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[1]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[0]*rdvSq4nuR; 
    outr[12] += 0.9682458365518543*vtSq[1]*fr[20]*rdvSq4nuR+0.3872983346207416*vtSq[1]*fr[16]*rdvSq4nuR-0.75*vtSq[3]*fr[13]*rdvSq4nuR-0.6708203932499369*vtSq[4]*fr[12]*rdvSq4nuR-0.75*vtSq[0]*fr[12]*rdvSq4nuR+0.4330127018922193*vtSq[2]*fr[6]*rdvSq4nuR-0.75*vtSq[1]*fr[5]*rdvSq4nuR+0.3872983346207416*fr[1]*vtSq[4]*rdvSq4nuR+0.4330127018922193*fr[2]*vtSq[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[1]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[1]*rdvSq4nuR; 
    outr[13] += 0.9682458365518543*vtSq[2]*fr[20]*rdvSq4nuR+0.3872983346207416*vtSq[2]*fr[17]*rdvSq4nuR-0.6708203932499369*vtSq[5]*fr[13]*rdvSq4nuR-0.75*vtSq[0]*fr[13]*rdvSq4nuR-0.75*vtSq[3]*fr[12]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[6]*rdvSq4nuR+0.3872983346207416*fr[2]*vtSq[5]*rdvSq4nuR-0.75*vtSq[2]*fr[5]*rdvSq4nuR+0.4330127018922193*fr[1]*vtSq[3]*rdvSq4nuR+0.4330127018922193*fr[0]*vtSq[2]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[2]*rdvSq4nuR; 
    outr[14] += (-0.75*vtSq[0]*fr[14]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[8]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[7]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[3]*rdvSq4nuR; 
    outr[15] += (-0.75*vtSq[0]*fr[15]*rdvSq4nuR)+0.4330127018922193*vtSq[2]*fr[10]*rdvSq4nuR+0.4330127018922193*vtSq[1]*fr[9]*rdvSq4nuR+0.4330127018922193*vtSq[0]*fr[4]*rdvSq4nuR; 
    outr[20] += (-3.75*vtSq[0]*fr[20]*rdvSq4nuR)-1.677050983124842*vtSq[5]*fr[17]*rdvSq4nuR-1.677050983124842*vtSq[4]*fr[16]*rdvSq4nuR+2.904737509655563*vtSq[2]*fr[13]*rdvSq4nuR+2.904737509655563*vtSq[1]*fr[12]*rdvSq4nuR-1.677050983124842*vtSq[3]*fr[6]*rdvSq4nuR+2.904737509655563*vtSq[0]*fr[5]*rdvSq4nuR-1.677050983124842*fr[2]*vtSq[2]*rdvSq4nuR-1.677050983124842*fr[1]*vtSq[1]*rdvSq4nuR-1.677050983124842*fr[0]*vtSq[0]*rdvSq4nuR; 

  } else {

    outl[5] += (-0.9682458365518543*vtSq[0]*fl[20]*rdvSq4nuL)-0.4330127018922193*vtSq[5]*fl[17]*rdvSq4nuL-0.4330127018922193*vtSq[4]*fl[16]*rdvSq4nuL-0.75*vtSq[2]*fl[13]*rdvSq4nuL-0.75*vtSq[1]*fl[12]*rdvSq4nuL-0.4330127018922193*vtSq[3]*fl[6]*rdvSq4nuL-0.75*vtSq[0]*fl[5]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[2]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[1]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[0]*rdvSq4nuL; 
    outl[12] += (-0.9682458365518543*vtSq[1]*fl[20]*rdvSq4nuL)-0.3872983346207416*vtSq[1]*fl[16]*rdvSq4nuL-0.75*vtSq[3]*fl[13]*rdvSq4nuL-0.6708203932499369*vtSq[4]*fl[12]*rdvSq4nuL-0.75*vtSq[0]*fl[12]*rdvSq4nuL-0.4330127018922193*vtSq[2]*fl[6]*rdvSq4nuL-0.75*vtSq[1]*fl[5]*rdvSq4nuL-0.3872983346207416*fl[1]*vtSq[4]*rdvSq4nuL-0.4330127018922193*fl[2]*vtSq[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[1]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[1]*rdvSq4nuL; 
    outl[13] += (-0.9682458365518543*vtSq[2]*fl[20]*rdvSq4nuL)-0.3872983346207416*vtSq[2]*fl[17]*rdvSq4nuL-0.6708203932499369*vtSq[5]*fl[13]*rdvSq4nuL-0.75*vtSq[0]*fl[13]*rdvSq4nuL-0.75*vtSq[3]*fl[12]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[6]*rdvSq4nuL-0.3872983346207416*fl[2]*vtSq[5]*rdvSq4nuL-0.75*vtSq[2]*fl[5]*rdvSq4nuL-0.4330127018922193*fl[1]*vtSq[3]*rdvSq4nuL-0.4330127018922193*fl[0]*vtSq[2]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[2]*rdvSq4nuL; 
    outl[14] += (-0.75*vtSq[0]*fl[14]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[8]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[7]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[3]*rdvSq4nuL; 
    outl[15] += (-0.75*vtSq[0]*fl[15]*rdvSq4nuL)-0.4330127018922193*vtSq[2]*fl[10]*rdvSq4nuL-0.4330127018922193*vtSq[1]*fl[9]*rdvSq4nuL-0.4330127018922193*vtSq[0]*fl[4]*rdvSq4nuL; 
    outl[20] += (-3.75*vtSq[0]*fl[20]*rdvSq4nuL)-1.677050983124842*vtSq[5]*fl[17]*rdvSq4nuL-1.677050983124842*vtSq[4]*fl[16]*rdvSq4nuL-2.904737509655563*vtSq[2]*fl[13]*rdvSq4nuL-2.904737509655563*vtSq[1]*fl[12]*rdvSq4nuL-1.677050983124842*vtSq[3]*fl[6]*rdvSq4nuL-2.904737509655563*vtSq[0]*fl[5]*rdvSq4nuL-1.677050983124842*fl[2]*vtSq[2]*rdvSq4nuL-1.677050983124842*fl[1]*vtSq[1]*rdvSq4nuL-1.677050983124842*fl[0]*vtSq[0]*rdvSq4nuL; 

  }
  return 0.0; 
} 
