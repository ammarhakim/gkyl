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

  double incr2[6]; 
  if (idxr[2] == 1) {

    incr2[3] = vtSq[0]*(0.4330127018922193*fr[0]-0.75*fr[3])+0.4330127018922193*(fr[2]*vtSq[2]+fr[1]*vtSq[1]); 

    outr[3] += incr2[3]*rdvSq4nuR; 

  } else {

    incr2[3] = vtSq[0]*((-0.75*fl[3])-0.4330127018922193*fl[0])-0.4330127018922193*(fl[2]*vtSq[2]+fl[1]*vtSq[1]); 

    outl[3] += incr2[3]*rdvSq4nuL; 

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

  double incr2[21]; 
  if (idxr[2] == 1) {

    incr2[3] = vtSq[0]*(0.9682458365518543*fr[18]-0.75*fr[3]+0.4330127018922193*fr[0])+0.4330127018922193*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(0.4330127018922193*fr[2]-0.75*fr[8])+vtSq[1]*(0.4330127018922193*fr[1]-0.75*fr[7])+0.4330127018922193*vtSq[3]*fr[6]; 
    incr2[7] = vtSq[1]*(0.9682458365518543*fr[18]+0.3872983346207416*fr[16]-0.75*fr[3]+0.4330127018922193*fr[0])+vtSq[3]*(0.4330127018922193*fr[2]-0.75*fr[8])+vtSq[4]*(0.3872983346207416*fr[1]-0.6708203932499369*fr[7])+vtSq[0]*(0.4330127018922193*fr[1]-0.75*fr[7])+0.4330127018922193*vtSq[2]*fr[6]; 
    incr2[8] = vtSq[2]*(0.9682458365518543*fr[18]+0.3872983346207416*fr[17]-0.75*fr[3]+0.4330127018922193*fr[0])+vtSq[5]*(0.3872983346207416*fr[2]-0.6708203932499369*fr[8])+vtSq[0]*(0.4330127018922193*fr[2]-0.75*fr[8])+vtSq[3]*(0.4330127018922193*fr[1]-0.75*fr[7])+0.4330127018922193*vtSq[1]*fr[6]; 
    incr2[11] = vtSq[0]*(0.4330127018922193*fr[4]-0.75*fr[11])+0.4330127018922193*(vtSq[2]*fr[10]+vtSq[1]*fr[9]); 
    incr2[14] = vtSq[0]*(0.4330127018922193*fr[5]-0.75*fr[14])+0.4330127018922193*(vtSq[2]*fr[13]+vtSq[1]*fr[12]); 
    incr2[18] = vtSq[0]*((-3.75*fr[18])+2.904737509655563*fr[3]-1.677050983124842*fr[0])-1.677050983124842*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(2.904737509655563*fr[8]-1.677050983124842*fr[2])+vtSq[1]*(2.904737509655563*fr[7]-1.677050983124842*fr[1])-1.677050983124842*vtSq[3]*fr[6]; 

    outr[3] += incr2[3]*rdvSq4nuR; 
    outr[7] += incr2[7]*rdvSq4nuR; 
    outr[8] += incr2[8]*rdvSq4nuR; 
    outr[11] += incr2[11]*rdvSq4nuR; 
    outr[14] += incr2[14]*rdvSq4nuR; 
    outr[18] += incr2[18]*rdvSq4nuR; 

  } else {

    incr2[3] = vtSq[0]*((-0.9682458365518543*fl[18])-0.75*fl[3]-0.4330127018922193*fl[0])-0.4330127018922193*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-0.75*fl[8])-0.4330127018922193*fl[2])+vtSq[1]*((-0.75*fl[7])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[3]*fl[6]; 
    incr2[7] = vtSq[1]*((-0.9682458365518543*fl[18])-0.3872983346207416*fl[16]-0.75*fl[3]-0.4330127018922193*fl[0])+vtSq[3]*((-0.75*fl[8])-0.4330127018922193*fl[2])+vtSq[4]*((-0.6708203932499369*fl[7])-0.3872983346207416*fl[1])+vtSq[0]*((-0.75*fl[7])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[2]*fl[6]; 
    incr2[8] = vtSq[2]*((-0.9682458365518543*fl[18])-0.3872983346207416*fl[17]-0.75*fl[3]-0.4330127018922193*fl[0])+vtSq[5]*((-0.6708203932499369*fl[8])-0.3872983346207416*fl[2])+vtSq[0]*((-0.75*fl[8])-0.4330127018922193*fl[2])+vtSq[3]*((-0.75*fl[7])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[1]*fl[6]; 
    incr2[11] = vtSq[0]*((-0.75*fl[11])-0.4330127018922193*fl[4])-0.4330127018922193*(vtSq[2]*fl[10]+vtSq[1]*fl[9]); 
    incr2[14] = vtSq[0]*((-0.75*fl[14])-0.4330127018922193*fl[5])-0.4330127018922193*(vtSq[2]*fl[13]+vtSq[1]*fl[12]); 
    incr2[18] = vtSq[0]*((-3.75*fl[18])-2.904737509655563*fl[3]-1.677050983124842*fl[0])-1.677050983124842*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-2.904737509655563*fl[8])-1.677050983124842*fl[2])+vtSq[1]*((-2.904737509655563*fl[7])-1.677050983124842*fl[1])-1.677050983124842*vtSq[3]*fl[6]; 

    outl[3] += incr2[3]*rdvSq4nuL; 
    outl[7] += incr2[7]*rdvSq4nuL; 
    outl[8] += incr2[8]*rdvSq4nuL; 
    outl[11] += incr2[11]*rdvSq4nuL; 
    outl[14] += incr2[14]*rdvSq4nuL; 
    outl[18] += incr2[18]*rdvSq4nuL; 

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

  double incr2[6]; 
  if (idxr[3] == 1) {

    incr2[4] = vtSq[0]*(0.4330127018922193*fr[0]-0.75*fr[4])+0.4330127018922193*(fr[2]*vtSq[2]+fr[1]*vtSq[1]); 

    outr[4] += incr2[4]*rdvSq4nuR; 

  } else {

    incr2[4] = vtSq[0]*((-0.75*fl[4])-0.4330127018922193*fl[0])-0.4330127018922193*(fl[2]*vtSq[2]+fl[1]*vtSq[1]); 

    outl[4] += incr2[4]*rdvSq4nuL; 

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

  double incr2[21]; 
  if (idxr[3] == 1) {

    incr2[4] = vtSq[0]*(0.9682458365518543*fr[19]-0.75*fr[4]+0.4330127018922193*fr[0])+0.4330127018922193*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(0.4330127018922193*fr[2]-0.75*fr[10])+vtSq[1]*(0.4330127018922193*fr[1]-0.75*fr[9])+0.4330127018922193*vtSq[3]*fr[6]; 
    incr2[9] = vtSq[1]*(0.9682458365518543*fr[19]+0.3872983346207416*fr[16]-0.75*fr[4]+0.4330127018922193*fr[0])+vtSq[3]*(0.4330127018922193*fr[2]-0.75*fr[10])+vtSq[4]*(0.3872983346207416*fr[1]-0.6708203932499369*fr[9])+vtSq[0]*(0.4330127018922193*fr[1]-0.75*fr[9])+0.4330127018922193*vtSq[2]*fr[6]; 
    incr2[10] = vtSq[2]*(0.9682458365518543*fr[19]+0.3872983346207416*fr[17]-0.75*fr[4]+0.4330127018922193*fr[0])+vtSq[5]*(0.3872983346207416*fr[2]-0.6708203932499369*fr[10])+vtSq[0]*(0.4330127018922193*fr[2]-0.75*fr[10])+vtSq[3]*(0.4330127018922193*fr[1]-0.75*fr[9])+0.4330127018922193*vtSq[1]*fr[6]; 
    incr2[11] = vtSq[0]*(0.4330127018922193*fr[3]-0.75*fr[11])+0.4330127018922193*(vtSq[2]*fr[8]+vtSq[1]*fr[7]); 
    incr2[15] = vtSq[0]*(0.4330127018922193*fr[5]-0.75*fr[15])+0.4330127018922193*(vtSq[2]*fr[13]+vtSq[1]*fr[12]); 
    incr2[19] = vtSq[0]*((-3.75*fr[19])+2.904737509655563*fr[4]-1.677050983124842*fr[0])-1.677050983124842*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(2.904737509655563*fr[10]-1.677050983124842*fr[2])+vtSq[1]*(2.904737509655563*fr[9]-1.677050983124842*fr[1])-1.677050983124842*vtSq[3]*fr[6]; 

    outr[4] += incr2[4]*rdvSq4nuR; 
    outr[9] += incr2[9]*rdvSq4nuR; 
    outr[10] += incr2[10]*rdvSq4nuR; 
    outr[11] += incr2[11]*rdvSq4nuR; 
    outr[15] += incr2[15]*rdvSq4nuR; 
    outr[19] += incr2[19]*rdvSq4nuR; 

  } else {

    incr2[4] = vtSq[0]*((-0.9682458365518543*fl[19])-0.75*fl[4]-0.4330127018922193*fl[0])-0.4330127018922193*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-0.75*fl[10])-0.4330127018922193*fl[2])+vtSq[1]*((-0.75*fl[9])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[3]*fl[6]; 
    incr2[9] = vtSq[1]*((-0.9682458365518543*fl[19])-0.3872983346207416*fl[16]-0.75*fl[4]-0.4330127018922193*fl[0])+vtSq[3]*((-0.75*fl[10])-0.4330127018922193*fl[2])+vtSq[4]*((-0.6708203932499369*fl[9])-0.3872983346207416*fl[1])+vtSq[0]*((-0.75*fl[9])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[2]*fl[6]; 
    incr2[10] = vtSq[2]*((-0.9682458365518543*fl[19])-0.3872983346207416*fl[17]-0.75*fl[4]-0.4330127018922193*fl[0])+vtSq[5]*((-0.6708203932499369*fl[10])-0.3872983346207416*fl[2])+vtSq[0]*((-0.75*fl[10])-0.4330127018922193*fl[2])+vtSq[3]*((-0.75*fl[9])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[1]*fl[6]; 
    incr2[11] = vtSq[0]*((-0.75*fl[11])-0.4330127018922193*fl[3])-0.4330127018922193*(vtSq[2]*fl[8]+vtSq[1]*fl[7]); 
    incr2[15] = vtSq[0]*((-0.75*fl[15])-0.4330127018922193*fl[5])-0.4330127018922193*(vtSq[2]*fl[13]+vtSq[1]*fl[12]); 
    incr2[19] = vtSq[0]*((-3.75*fl[19])-2.904737509655563*fl[4]-1.677050983124842*fl[0])-1.677050983124842*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-2.904737509655563*fl[10])-1.677050983124842*fl[2])+vtSq[1]*((-2.904737509655563*fl[9])-1.677050983124842*fl[1])-1.677050983124842*vtSq[3]*fl[6]; 

    outl[4] += incr2[4]*rdvSq4nuL; 
    outl[9] += incr2[9]*rdvSq4nuL; 
    outl[10] += incr2[10]*rdvSq4nuL; 
    outl[11] += incr2[11]*rdvSq4nuL; 
    outl[15] += incr2[15]*rdvSq4nuL; 
    outl[19] += incr2[19]*rdvSq4nuL; 

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

  double incr2[6]; 
  if (idxr[4] == 1) {

    incr2[5] = vtSq[0]*(0.4330127018922193*fr[0]-0.75*fr[5])+0.4330127018922193*(fr[2]*vtSq[2]+fr[1]*vtSq[1]); 

    outr[5] += incr2[5]*rdvSq4nuR; 

  } else {

    incr2[5] = vtSq[0]*((-0.75*fl[5])-0.4330127018922193*fl[0])-0.4330127018922193*(fl[2]*vtSq[2]+fl[1]*vtSq[1]); 

    outl[5] += incr2[5]*rdvSq4nuL; 

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

  double incr2[21]; 
  if (idxr[4] == 1) {

    incr2[5] = vtSq[0]*(0.9682458365518543*fr[20]-0.75*fr[5]+0.4330127018922193*fr[0])+0.4330127018922193*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(0.4330127018922193*fr[2]-0.75*fr[13])+vtSq[1]*(0.4330127018922193*fr[1]-0.75*fr[12])+0.4330127018922193*vtSq[3]*fr[6]; 
    incr2[12] = vtSq[1]*(0.9682458365518543*fr[20]+0.3872983346207416*fr[16]-0.75*fr[5]+0.4330127018922193*fr[0])+vtSq[3]*(0.4330127018922193*fr[2]-0.75*fr[13])+vtSq[4]*(0.3872983346207416*fr[1]-0.6708203932499369*fr[12])+vtSq[0]*(0.4330127018922193*fr[1]-0.75*fr[12])+0.4330127018922193*vtSq[2]*fr[6]; 
    incr2[13] = vtSq[2]*(0.9682458365518543*fr[20]+0.3872983346207416*fr[17]-0.75*fr[5]+0.4330127018922193*fr[0])+vtSq[5]*(0.3872983346207416*fr[2]-0.6708203932499369*fr[13])+vtSq[0]*(0.4330127018922193*fr[2]-0.75*fr[13])+vtSq[3]*(0.4330127018922193*fr[1]-0.75*fr[12])+0.4330127018922193*vtSq[1]*fr[6]; 
    incr2[14] = vtSq[0]*(0.4330127018922193*fr[3]-0.75*fr[14])+0.4330127018922193*(vtSq[2]*fr[8]+vtSq[1]*fr[7]); 
    incr2[15] = vtSq[0]*(0.4330127018922193*fr[4]-0.75*fr[15])+0.4330127018922193*(vtSq[2]*fr[10]+vtSq[1]*fr[9]); 
    incr2[20] = vtSq[0]*((-3.75*fr[20])+2.904737509655563*fr[5]-1.677050983124842*fr[0])-1.677050983124842*(vtSq[5]*fr[17]+vtSq[4]*fr[16])+vtSq[2]*(2.904737509655563*fr[13]-1.677050983124842*fr[2])+vtSq[1]*(2.904737509655563*fr[12]-1.677050983124842*fr[1])-1.677050983124842*vtSq[3]*fr[6]; 

    outr[5] += incr2[5]*rdvSq4nuR; 
    outr[12] += incr2[12]*rdvSq4nuR; 
    outr[13] += incr2[13]*rdvSq4nuR; 
    outr[14] += incr2[14]*rdvSq4nuR; 
    outr[15] += incr2[15]*rdvSq4nuR; 
    outr[20] += incr2[20]*rdvSq4nuR; 

  } else {

    incr2[5] = vtSq[0]*((-0.9682458365518543*fl[20])-0.75*fl[5]-0.4330127018922193*fl[0])-0.4330127018922193*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-0.75*fl[13])-0.4330127018922193*fl[2])+vtSq[1]*((-0.75*fl[12])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[3]*fl[6]; 
    incr2[12] = vtSq[1]*((-0.9682458365518543*fl[20])-0.3872983346207416*fl[16]-0.75*fl[5]-0.4330127018922193*fl[0])+vtSq[3]*((-0.75*fl[13])-0.4330127018922193*fl[2])+vtSq[4]*((-0.6708203932499369*fl[12])-0.3872983346207416*fl[1])+vtSq[0]*((-0.75*fl[12])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[2]*fl[6]; 
    incr2[13] = vtSq[2]*((-0.9682458365518543*fl[20])-0.3872983346207416*fl[17]-0.75*fl[5]-0.4330127018922193*fl[0])+vtSq[5]*((-0.6708203932499369*fl[13])-0.3872983346207416*fl[2])+vtSq[0]*((-0.75*fl[13])-0.4330127018922193*fl[2])+vtSq[3]*((-0.75*fl[12])-0.4330127018922193*fl[1])-0.4330127018922193*vtSq[1]*fl[6]; 
    incr2[14] = vtSq[0]*((-0.75*fl[14])-0.4330127018922193*fl[3])-0.4330127018922193*(vtSq[2]*fl[8]+vtSq[1]*fl[7]); 
    incr2[15] = vtSq[0]*((-0.75*fl[15])-0.4330127018922193*fl[4])-0.4330127018922193*(vtSq[2]*fl[10]+vtSq[1]*fl[9]); 
    incr2[20] = vtSq[0]*((-3.75*fl[20])-2.904737509655563*fl[5]-1.677050983124842*fl[0])-1.677050983124842*(vtSq[5]*fl[17]+vtSq[4]*fl[16])+vtSq[2]*((-2.904737509655563*fl[13])-1.677050983124842*fl[2])+vtSq[1]*((-2.904737509655563*fl[12])-1.677050983124842*fl[1])-1.677050983124842*vtSq[3]*fl[6]; 

    outl[5] += incr2[5]*rdvSq4nuL; 
    outl[12] += incr2[12]*rdvSq4nuL; 
    outl[13] += incr2[13]*rdvSq4nuL; 
    outl[14] += incr2[14]*rdvSq4nuL; 
    outl[15] += incr2[15]*rdvSq4nuL; 
    outl[20] += incr2[20]*rdvSq4nuL; 

  }
  return 0.0; 
} 
