#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (1.732050807568877*fr[0]-3.0*fr[1])*rdxSq2nur; 
    outr[4] += (1.732050807568877*fr[2]-3.0*fr[4])*rdxSq2nur; 
    outr[5] += (1.732050807568877*fr[3]-3.0*fr[5])*rdxSq2nur; 
    outr[7] += (1.732050807568877*fr[6]-3.0*fr[7])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.0*fl[1])-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[4] += ((-3.0*fl[4])-1.732050807568877*fl[2])*rdxSq2nul; 
    outl[5] += ((-3.0*fl[5])-1.732050807568877*fl[3])*rdxSq2nul; 
    outl[7] += ((-3.0*fl[7])-1.732050807568877*fl[6])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (3.872983346207417*fr[7]-3.0*fr[1]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[4] += (3.872983346207417*fr[11]-3.0*fr[4]+1.732050807568877*fr[2])*rdxSq2nur; 
    outr[5] += (3.872983346207417*fr[13]-3.0*fr[5]+1.732050807568877*fr[3])*rdxSq2nur; 
    outr[7] += ((-15.0*fr[7])+11.61895003862225*fr[1]-6.708203932499369*fr[0])*rdxSq2nur; 
    outr[10] += (3.872983346207417*fr[17]-3.0*fr[10]+1.732050807568877*fr[6])*rdxSq2nur; 
    outr[11] += ((-15.0*fr[11])+11.61895003862225*fr[4]-6.708203932499369*fr[2])*rdxSq2nur; 
    outr[12] += (1.732050807568877*fr[8]-3.0*fr[12])*rdxSq2nur; 
    outr[13] += ((-15.0*fr[13])+11.61895003862225*fr[5]-6.708203932499369*fr[3])*rdxSq2nur; 
    outr[15] += (1.732050807568877*fr[9]-3.0*fr[15])*rdxSq2nur; 
    outr[17] += ((-15.0*fr[17])+11.61895003862225*fr[10]-6.708203932499369*fr[6])*rdxSq2nur; 
    outr[18] += (1.732050807568877*fr[14]-3.0*fr[18])*rdxSq2nur; 
    outr[19] += (1.732050807568877*fr[16]-3.0*fr[19])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.872983346207417*fl[7])-3.0*fl[1]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[4] += ((-3.872983346207417*fl[11])-3.0*fl[4]-1.732050807568877*fl[2])*rdxSq2nul; 
    outl[5] += ((-3.872983346207417*fl[13])-3.0*fl[5]-1.732050807568877*fl[3])*rdxSq2nul; 
    outl[7] += ((-15.0*fl[7])-11.61895003862225*fl[1]-6.708203932499369*fl[0])*rdxSq2nul; 
    outl[10] += ((-3.872983346207417*fl[17])-3.0*fl[10]-1.732050807568877*fl[6])*rdxSq2nul; 
    outl[11] += ((-15.0*fl[11])-11.61895003862225*fl[4]-6.708203932499369*fl[2])*rdxSq2nul; 
    outl[12] += ((-3.0*fl[12])-1.732050807568877*fl[8])*rdxSq2nul; 
    outl[13] += ((-15.0*fl[13])-11.61895003862225*fl[5]-6.708203932499369*fl[3])*rdxSq2nul; 
    outl[15] += ((-3.0*fl[15])-1.732050807568877*fl[9])*rdxSq2nul; 
    outl[17] += ((-15.0*fl[17])-11.61895003862225*fl[10]-6.708203932499369*fl[6])*rdxSq2nul; 
    outl[18] += ((-3.0*fl[18])-1.732050807568877*fl[14])*rdxSq2nul; 
    outl[19] += ((-3.0*fl[19])-1.732050807568877*fl[16])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  if (idxr[2] == 1) {

    outr[2] += (1.732050807568877*fr[0]-3.0*fr[2])*rdxSq2nur; 
    outr[4] += (1.732050807568877*fr[1]-3.0*fr[4])*rdxSq2nur; 
    outr[6] += (1.732050807568877*fr[3]-3.0*fr[6])*rdxSq2nur; 
    outr[7] += (1.732050807568877*fr[5]-3.0*fr[7])*rdxSq2nur; 

  } else {

    outl[2] += ((-3.0*fl[2])-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[4] += ((-3.0*fl[4])-1.732050807568877*fl[1])*rdxSq2nul; 
    outl[6] += ((-3.0*fl[6])-1.732050807568877*fl[3])*rdxSq2nul; 
    outl[7] += ((-3.0*fl[7])-1.732050807568877*fl[5])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  if (idxr[2] == 1) {

    outr[2] += (3.872983346207417*fr[8]-3.0*fr[2]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[4] += (3.872983346207417*fr[12]-3.0*fr[4]+1.732050807568877*fr[1])*rdxSq2nur; 
    outr[6] += (3.872983346207417*fr[14]-3.0*fr[6]+1.732050807568877*fr[3])*rdxSq2nur; 
    outr[8] += ((-15.0*fr[8])+11.61895003862225*fr[2]-6.708203932499369*fr[0])*rdxSq2nur; 
    outr[10] += (3.872983346207417*fr[18]-3.0*fr[10]+1.732050807568877*fr[5])*rdxSq2nur; 
    outr[11] += (1.732050807568877*fr[7]-3.0*fr[11])*rdxSq2nur; 
    outr[12] += ((-15.0*fr[12])+11.61895003862225*fr[4]-6.708203932499369*fr[1])*rdxSq2nur; 
    outr[14] += ((-15.0*fr[14])+11.61895003862225*fr[6]-6.708203932499369*fr[3])*rdxSq2nur; 
    outr[16] += (1.732050807568877*fr[9]-3.0*fr[16])*rdxSq2nur; 
    outr[17] += (1.732050807568877*fr[13]-3.0*fr[17])*rdxSq2nur; 
    outr[18] += ((-15.0*fr[18])+11.61895003862225*fr[10]-6.708203932499369*fr[5])*rdxSq2nur; 
    outr[19] += (1.732050807568877*fr[15]-3.0*fr[19])*rdxSq2nur; 

  } else {

    outl[2] += ((-3.872983346207417*fl[8])-3.0*fl[2]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[4] += ((-3.872983346207417*fl[12])-3.0*fl[4]-1.732050807568877*fl[1])*rdxSq2nul; 
    outl[6] += ((-3.872983346207417*fl[14])-3.0*fl[6]-1.732050807568877*fl[3])*rdxSq2nul; 
    outl[8] += ((-15.0*fl[8])-11.61895003862225*fl[2]-6.708203932499369*fl[0])*rdxSq2nul; 
    outl[10] += ((-3.872983346207417*fl[18])-3.0*fl[10]-1.732050807568877*fl[5])*rdxSq2nul; 
    outl[11] += ((-3.0*fl[11])-1.732050807568877*fl[7])*rdxSq2nul; 
    outl[12] += ((-15.0*fl[12])-11.61895003862225*fl[4]-6.708203932499369*fl[1])*rdxSq2nul; 
    outl[14] += ((-15.0*fl[14])-11.61895003862225*fl[6]-6.708203932499369*fl[3])*rdxSq2nul; 
    outl[16] += ((-3.0*fl[16])-1.732050807568877*fl[9])*rdxSq2nul; 
    outl[17] += ((-3.0*fl[17])-1.732050807568877*fl[13])*rdxSq2nul; 
    outl[18] += ((-15.0*fl[18])-11.61895003862225*fl[10]-6.708203932499369*fl[5])*rdxSq2nul; 
    outl[19] += ((-3.0*fl[19])-1.732050807568877*fl[15])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 

  if (idxr[3] == 1) {

    outr[3] += (1.732050807568877*fr[0]-3.0*fr[3])*rdxSq2nur; 
    outr[5] += (1.732050807568877*fr[1]-3.0*fr[5])*rdxSq2nur; 
    outr[6] += (1.732050807568877*fr[2]-3.0*fr[6])*rdxSq2nur; 
    outr[7] += (1.732050807568877*fr[4]-3.0*fr[7])*rdxSq2nur; 

  } else {

    outl[3] += ((-3.0*fl[3])-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[5] += ((-3.0*fl[5])-1.732050807568877*fl[1])*rdxSq2nul; 
    outl[6] += ((-3.0*fl[6])-1.732050807568877*fl[2])*rdxSq2nul; 
    outl[7] += ((-3.0*fl[7])-1.732050807568877*fl[4])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq2nur = 2.0*nu[2]/(dxvr[2]*dxvr[2]); 

  if (idxr[3] == 1) {

    outr[3] += (3.872983346207417*fr[9]-3.0*fr[3]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[5] += (3.872983346207417*fr[15]-3.0*fr[5]+1.732050807568877*fr[1])*rdxSq2nur; 
    outr[6] += (3.872983346207417*fr[16]-3.0*fr[6]+1.732050807568877*fr[2])*rdxSq2nur; 
    outr[9] += ((-15.0*fr[9])+11.61895003862225*fr[3]-6.708203932499369*fr[0])*rdxSq2nur; 
    outr[10] += (3.872983346207417*fr[19]-3.0*fr[10]+1.732050807568877*fr[4])*rdxSq2nur; 
    outr[13] += (1.732050807568877*fr[7]-3.0*fr[13])*rdxSq2nur; 
    outr[14] += (1.732050807568877*fr[8]-3.0*fr[14])*rdxSq2nur; 
    outr[15] += ((-15.0*fr[15])+11.61895003862225*fr[5]-6.708203932499369*fr[1])*rdxSq2nur; 
    outr[16] += ((-15.0*fr[16])+11.61895003862225*fr[6]-6.708203932499369*fr[2])*rdxSq2nur; 
    outr[17] += (1.732050807568877*fr[11]-3.0*fr[17])*rdxSq2nur; 
    outr[18] += (1.732050807568877*fr[12]-3.0*fr[18])*rdxSq2nur; 
    outr[19] += ((-15.0*fr[19])+11.61895003862225*fr[10]-6.708203932499369*fr[4])*rdxSq2nur; 

  } else {

    outl[3] += ((-3.872983346207417*fl[9])-3.0*fl[3]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[5] += ((-3.872983346207417*fl[15])-3.0*fl[5]-1.732050807568877*fl[1])*rdxSq2nul; 
    outl[6] += ((-3.872983346207417*fl[16])-3.0*fl[6]-1.732050807568877*fl[2])*rdxSq2nul; 
    outl[9] += ((-15.0*fl[9])-11.61895003862225*fl[3]-6.708203932499369*fl[0])*rdxSq2nul; 
    outl[10] += ((-3.872983346207417*fl[19])-3.0*fl[10]-1.732050807568877*fl[4])*rdxSq2nul; 
    outl[13] += ((-3.0*fl[13])-1.732050807568877*fl[7])*rdxSq2nul; 
    outl[14] += ((-3.0*fl[14])-1.732050807568877*fl[8])*rdxSq2nul; 
    outl[15] += ((-15.0*fl[15])-11.61895003862225*fl[5]-6.708203932499369*fl[1])*rdxSq2nul; 
    outl[16] += ((-15.0*fl[16])-11.61895003862225*fl[6]-6.708203932499369*fl[2])*rdxSq2nul; 
    outl[17] += ((-3.0*fl[17])-1.732050807568877*fl[11])*rdxSq2nul; 
    outl[18] += ((-3.0*fl[18])-1.732050807568877*fl[12])*rdxSq2nul; 
    outl[19] += ((-15.0*fl[19])-11.61895003862225*fl[10]-6.708203932499369*fl[4])*rdxSq2nul; 

  }

} 
