#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xMaxP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 

} 
