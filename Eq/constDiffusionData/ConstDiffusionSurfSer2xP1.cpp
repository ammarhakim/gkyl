#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf2xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 

} 
void ConstDiffusionSurf2xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf2xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  incr1[0] = (-1.623797632095822*fr[1])-1.623797632095822*fl[1]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = 2.8125*fr[1]+2.8125*fl[1]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[2] = (-1.623797632095822*fr[3])-1.623797632095822*fl[3]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = 2.8125*fr[3]+2.8125*fl[3]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 

  incr2[1] = 0.75*fr[1]-0.75*fl[1]; 
  incr2[3] = 0.75*fr[3]-0.75*fl[3]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr2[3]*rdxFnur)-1.0*incr1[3]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf2xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  incr1[0] = (-1.623797632095822*fr[2])-1.623797632095822*fl[2]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[3])-1.623797632095822*fl[3]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = 2.8125*fr[2]+2.8125*fl[2]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[3] = 2.8125*fr[3]+2.8125*fl[3]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 

  incr2[2] = 0.75*fr[2]-0.75*fl[2]; 
  incr2[3] = 0.75*fr[3]-0.75*fl[3]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr2[2]*rdxFnur)-1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr2[3]*rdxFnur)-1.0*incr1[3]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 

} 
