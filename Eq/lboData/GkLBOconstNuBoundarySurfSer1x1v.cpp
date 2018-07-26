#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[1*2]: bulk velocity (in 1 directions) times nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*nuVtSq[0]; 
  alpha[1] = 1.414213562373095*nuVtSq[1]; 

  if (idxr[1] == 1) {

    outr[2] += (0.4330127018922193*(alpha[1]*fr[1]+alpha[0]*fr[0])-0.75*(alpha[1]*fr[3]+alpha[0]*fr[2]))*rdvSq4r; 
    outr[3] += (0.4330127018922193*(alpha[0]*fr[1]+fr[0]*alpha[1])-0.75*(alpha[0]*fr[3]+alpha[1]*fr[2]))*rdvSq4r; 

  } else {

    outl[2] += ((-0.75*(alpha[1]*fl[3]+alpha[0]*fl[2]))-0.4330127018922193*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4l; 
    outl[3] += ((-0.75*(alpha[0]*fl[3]+alpha[1]*fl[2]))-0.4330127018922193*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4l; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P2(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[1*3]: bulk velocity (in 1 directions) times nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*nuVtSq[0]; 
  alpha[1] = 1.414213562373095*nuVtSq[1]; 
  alpha[4] = 1.414213562373095*nuVtSq[2]; 

  if (idxr[1] == 1) {

    outr[2] += (0.9682458365518543*alpha[1]*fr[7]-0.75*alpha[4]*fr[6]+0.9682458365518543*alpha[0]*fr[5]+0.4330127018922193*alpha[4]*fr[4]-0.75*(alpha[1]*fr[3]+alpha[0]*fr[2])+0.4330127018922193*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4r; 
    outr[3] += ((0.8660254037844387*alpha[4]+0.9682458365518543*alpha[0])*fr[7]+alpha[1]*((-0.6708203932499369*fr[6])+0.9682458365518543*fr[5]+0.3872983346207416*fr[4])+(0.3872983346207416*fr[1]-0.6708203932499369*fr[3])*alpha[4]-0.75*(alpha[0]*fr[3]+alpha[1]*fr[2])+0.4330127018922193*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4r; 
    outr[5] += ((-3.75*alpha[1]*fr[7])+2.904737509655563*alpha[4]*fr[6]-3.75*alpha[0]*fr[5]-1.677050983124842*alpha[4]*fr[4]+2.904737509655563*(alpha[1]*fr[3]+alpha[0]*fr[2])-1.677050983124842*(alpha[1]*fr[1]+alpha[0]*fr[0]))*rdvSq4r; 
    outr[6] += (0.8660254037844386*alpha[1]*fr[7]+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fr[6]+0.9682458365518543*alpha[4]*fr[5]+(0.276641667586244*alpha[4]+0.4330127018922194*alpha[0])*fr[4]+(0.4330127018922194*fr[0]-0.75*fr[2])*alpha[4]+alpha[1]*(0.3872983346207417*fr[1]-0.6708203932499369*fr[3]))*rdvSq4r; 
    outr[7] += (((-3.354101966249685*alpha[4])-3.75*alpha[0])*fr[7]+alpha[1]*(2.598076211353316*fr[6]-3.75*fr[5]-1.5*fr[4])+(2.598076211353316*fr[3]-1.5*fr[1])*alpha[4]+2.904737509655563*(alpha[0]*fr[3]+alpha[1]*fr[2])-1.677050983124842*(alpha[0]*fr[1]+fr[0]*alpha[1]))*rdvSq4r; 

  } else {

    outl[2] += ((-0.9682458365518543*alpha[1]*fl[7])-0.75*alpha[4]*fl[6]-0.9682458365518543*alpha[0]*fl[5]-0.4330127018922193*alpha[4]*fl[4]-0.75*(alpha[1]*fl[3]+alpha[0]*fl[2])-0.4330127018922193*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4l; 
    outl[3] += (((-0.8660254037844387*alpha[4])-0.9682458365518543*alpha[0])*fl[7]+alpha[1]*((-0.6708203932499369*fl[6])-0.9682458365518543*fl[5]-0.3872983346207416*fl[4])+((-0.6708203932499369*fl[3])-0.3872983346207416*fl[1])*alpha[4]-0.75*(alpha[0]*fl[3]+alpha[1]*fl[2])-0.4330127018922193*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4l; 
    outl[5] += ((-3.75*alpha[1]*fl[7])-2.904737509655563*alpha[4]*fl[6]-3.75*alpha[0]*fl[5]-1.677050983124842*alpha[4]*fl[4]-2.904737509655563*(alpha[1]*fl[3]+alpha[0]*fl[2])-1.677050983124842*(alpha[1]*fl[1]+alpha[0]*fl[0]))*rdvSq4l; 
    outl[6] += ((-0.8660254037844386*alpha[1]*fl[7])+((-0.479157423749955*alpha[4])-0.75*alpha[0])*fl[6]-0.9682458365518543*alpha[4]*fl[5]+((-0.276641667586244*alpha[4])-0.4330127018922194*alpha[0])*fl[4]+((-0.75*fl[2])-0.4330127018922194*fl[0])*alpha[4]+alpha[1]*((-0.6708203932499369*fl[3])-0.3872983346207417*fl[1]))*rdvSq4l; 
    outl[7] += (((-3.354101966249685*alpha[4])-3.75*alpha[0])*fl[7]+alpha[1]*((-2.598076211353316*fl[6])-3.75*fl[5]-1.5*fl[4])+((-2.598076211353316*fl[3])-1.5*fl[1])*alpha[4]-2.904737509655563*(alpha[0]*fl[3]+alpha[1]*fl[2])-1.677050983124842*(alpha[0]*fl[1]+fl[0]*alpha[1]))*rdvSq4l; 

  }
  return 0.0; 
} 
