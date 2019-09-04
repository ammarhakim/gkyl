#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dxv[1]:    Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (1.732050807568877*fr[0]-3.0*fr[1])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.0*fl[1])-1.732050807568877*fl[0])*rdxSq2nul; 

  }

} 
void ConstDiffusionBoundarySurf1xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dxv[1]:    Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  if (idxr[1] == 1) {

    outr[1] += (3.872983346207417*fr[2]-3.0*fr[1]+1.732050807568877*fr[0])*rdxSq2nur; 
    outr[2] += ((-15.0*fr[2])+11.61895003862225*fr[1]-6.708203932499369*fr[0])*rdxSq2nur; 

  } else {

    outl[1] += ((-3.872983346207417*fl[2])-3.0*fl[1]-1.732050807568877*fl[0])*rdxSq2nul; 
    outl[2] += ((-15.0*fl[2])-11.61895003862225*fl[1]-6.708203932499369*fl[0])*rdxSq2nul; 

  }

} 
