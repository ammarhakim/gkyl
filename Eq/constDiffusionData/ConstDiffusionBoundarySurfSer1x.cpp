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

  double incr1[2]; 
  double incr2[2]; 

  incr1[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr1[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  incr2[1] = (-1.0*fr[1])+fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 

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

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-1.341640786499874*fr[2])+1.341640786499874*fl[2]+2.381569860407206*fr[1]+2.381569860407206*fl[1]-1.875*fr[0]+1.875*fl[0]; 
  incr1[1] = 2.32379000772445*fr[2]-2.32379000772445*fl[2]-4.125*fr[1]-4.125*fl[1]+3.247595264191645*fr[0]-3.247595264191645*fl[0]; 
  incr1[2] = (-3.0*fr[2])+3.0*fl[2]+5.325352101035199*fr[1]+5.325352101035199*fl[1]-4.192627457812106*fr[0]+4.192627457812106*fl[0]; 

  incr2[1] = 0.8472151069828725*fr[2]+0.8472151069828725*fl[2]-1.21875*fr[1]+1.21875*fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
  incr2[2] = (-3.28125*fr[2])-3.28125*fl[2]+4.720198453190289*fr[1]-4.720198453190289*fl[1]-3.354101966249685*fr[0]-3.354101966249685*fl[0]; 

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
  outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
  outl[2] += incr2[2]*rdxSq2nul-1.0*incr1[2]*rdxSq2nul; 

} 
