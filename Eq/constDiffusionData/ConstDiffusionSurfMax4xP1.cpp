#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf4xMaxP1_diffDirs2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs23_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs23_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs234_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs234_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs234_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs1234_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs1234_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs1234_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs1234_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs123_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs123_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs123_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs24_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs24_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs124_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs124_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs124_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs12_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs12_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs34_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs34_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs134_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs134_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs134_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs13_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs13_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs4_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs14_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs14_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5625*fl[1]-0.5625*fr[1]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xMaxP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[5]; 
  double incr2[5]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5625*fl[2]-0.5625*fr[2]; 
  incr1[3] = 0.5625*fl[3]-0.5625*fr[3]; 
  incr1[4] = 0.5625*fl[4]-0.5625*fr[4]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 

} 
