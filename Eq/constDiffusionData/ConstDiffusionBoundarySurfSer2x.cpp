#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (1.732050807568877*fr[0]-3.0*fr[1])*rdxSq2nur; 
    outr[3] += (1.732050807568877*fr[2]-3.0*fr[3])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.0*fl[1])-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[3] += ((-3.0*fl[3])-1.732050807568877*fl[2])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf2xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (3.872983346207417*fr[4]-3.0*fr[1]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[3] += (3.872983346207417*fr[6]-3.0*fr[3]+1.732050807568877*fr[2])*rdxSq2nur; 
    outr[4] += ((-15.0*fr[4])+11.61895003862225*fr[1]-6.708203932499369*fr[0])*rdxSq2nur; 
    outr[6] += ((-15.0*fr[6])+11.61895003862225*fr[3]-6.708203932499369*fr[2])*rdxSq2nur; 
    outr[7] += (1.732050807568877*fr[5]-3.0*fr[7])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.872983346207417*fl[4])-3.0*fl[1]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[3] += ((-3.872983346207417*fl[6])-3.0*fl[3]-1.732050807568877*fl[2])*rdxSq2nul; 
    outl[4] += ((-15.0*fl[4])-11.61895003862225*fl[1]-6.708203932499369*fl[0])*rdxSq2nul; 
    outl[6] += ((-15.0*fl[6])-11.61895003862225*fl[3]-6.708203932499369*fl[2])*rdxSq2nul; 
    outl[7] += ((-3.0*fl[7])-1.732050807568877*fl[5])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  if (idxr[2] == 1) {

    outr[2] += (1.732050807568877*fr[0]-3.0*fr[2])*rdxSq2nur; 
    outr[3] += (1.732050807568877*fr[1]-3.0*fr[3])*rdxSq2nur; 

  } else {

    outl[2] += ((-3.0*fl[2])-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[3] += ((-3.0*fl[3])-1.732050807568877*fl[1])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf2xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  if (idxr[2] == 1) {

    outr[2] += (3.872983346207417*fr[5]-3.0*fr[2]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[3] += (3.872983346207417*fr[7]-3.0*fr[3]+1.732050807568877*fr[1])*rdxSq2nur; 
    outr[5] += ((-15.0*fr[5])+11.61895003862225*fr[2]-6.708203932499369*fr[0])*rdxSq2nur; 
    outr[6] += (1.732050807568877*fr[4]-3.0*fr[6])*rdxSq2nur; 
    outr[7] += ((-15.0*fr[7])+11.61895003862225*fr[3]-6.708203932499369*fr[1])*rdxSq2nur; 

  } else {

    outl[2] += ((-3.872983346207417*fl[5])-3.0*fl[2]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[3] += ((-3.872983346207417*fl[7])-3.0*fl[3]-1.732050807568877*fl[1])*rdxSq2nul; 
    outl[5] += ((-15.0*fl[5])-11.61895003862225*fl[2]-6.708203932499369*fl[0])*rdxSq2nul; 
    outl[6] += ((-3.0*fl[6])-1.732050807568877*fl[4])*rdxSq2nul; 
    outl[7] += ((-15.0*fl[7])-11.61895003862225*fl[3]-6.708203932499369*fl[1])*rdxSq2nul; 

  }

} 
