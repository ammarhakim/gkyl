#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xMaxP2_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-0.6708203932499369*fr[2])+0.6708203932499369*fl[2]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[2]-1.161895003862225*fl[2]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = (-1.5*fr[2])+1.5*fl[2]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 

  incr2[1] = 0.4236075534914363*fr[2]+0.4236075534914363*fl[2]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[2] = (-1.640625*fr[2])-1.640625*fl[2]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += incr2[2]*rdxSq4nul-1.0*incr1[2]*rdxSq4nul; 

} 
