#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf2xMaxP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 

} 
void ConstDiffusionSurf2xMaxP1_diffDirs12_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 

} 
void ConstDiffusionSurf2xMaxP1_diffDirs12_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 

} 
void ConstDiffusionSurf2xMaxP1_diffDirs2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 

} 
