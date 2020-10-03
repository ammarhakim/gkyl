#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf3xSerP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[6] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs12_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[6] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs12_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[5] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs123_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[6] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs123_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[5] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs123_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs13_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[6] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs13_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[5] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs23_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[5] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs23_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSerP1_diffDirs3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 

} 
