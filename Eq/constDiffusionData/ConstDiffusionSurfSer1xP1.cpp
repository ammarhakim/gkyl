#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 
  double incr3[2]; 
  double incr4[2]; 

  incr1[0] = (-1.623797632095822*fr[1])-1.623797632095822*fl[1]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = 2.8125*fr[1]+2.8125*fl[1]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 

  incr2[1] = 0.75*fr[1]-0.75*fl[1]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 

} 
void ConstDiffusionVarCoeffSurf1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[2]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  incr1[0] = 0.6629126073623879*fr[1]*nul[1]+0.6629126073623879*fl[1]*nul[1]-0.6889189901577683*fr[0]*nul[1]+0.6889189901577683*fl[0]*nul[1]+0.3827327723098713*nul[0]*fr[1]+0.3827327723098713*nul[0]*fl[1]-0.3977475644174328*fr[0]*nul[0]+0.3977475644174328*fl[0]*nul[0]; 
  incr1[1] = (-1.148198316929614*fr[1]*nul[1])-1.148198316929614*fl[1]*nul[1]+1.193242693252298*fr[0]*nul[1]-1.193242693252298*fl[0]*nul[1]-0.6629126073623879*nul[0]*fr[1]-0.6629126073623879*nul[0]*fl[1]+0.6889189901577683*fr[0]*nul[0]-0.6889189901577683*fl[0]*nul[0]; 

  incr2[1] = (-0.6123724356957944*fr[1]*nul[1])+0.6123724356957944*fl[1]*nul[1]+0.5303300858899105*fr[0]*nul[1]+0.5303300858899105*fl[0]*nul[1]-0.3535533905932737*nul[0]*fr[1]+0.3535533905932737*nul[0]*fl[1]+0.3061862178478971*fr[0]*nul[0]+0.3061862178478971*fl[0]*nul[0]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr2[1]*rdxFr+incr1[1]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += incr1[1]*rdxFl-1.0*incr2[1]*rdxFl; 

} 
