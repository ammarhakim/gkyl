#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xMaxP2_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-0.6708203932499369*fr[2])+0.6708203932499369*fl[2]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[2]-1.161895003862225*fl[2]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = (-1.5*fr[2])+1.5*fl[2]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 

  incr2[1] = 0.4236075534914363*fr[2]+0.4236075534914363*fl[2]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[2] = (-1.640625*fr[2])-1.640625*fl[2]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf1xMaxP2_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = 6.708203932499369*fr[2]-6.708203932499369*fl[2]-8.11898816047911*fr[1]-8.11898816047911*fl[1]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-11.61895003862225*fr[2])+11.61895003862225*fl[2]+14.0625*fr[1]+14.0625*fl[1]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[2] = 10.5*fr[2]-10.5*fl[2]-10.16658128379447*fr[1]-10.16658128379447*fl[1]+4.192627457812107*fr[0]-4.192627457812107*fl[0]; 

  incr2[1] = (-2.541645320948617*fr[2])-2.541645320948617*fl[2]+1.40625*fr[1]-1.40625*fl[1]; 
  incr2[2] = 9.84375*fr[2]+9.84375*fl[2]-5.446382830604179*fr[1]+5.446382830604179*fl[1]; 

  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr2[2]*rdxFnur)-1.0*incr1[2]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf1xMaxP2_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  incr1[0] = (-35.21807064562169*fr[2])+35.21807064562169*fl[2]+34.09975027401226*fr[1]+34.09975027401226*fl[1]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 60.9994877027668*fr[2]-60.9994877027668*fl[2]-59.0625*fr[1]-59.0625*fl[1]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[2] = (-33.75*fr[2])+33.75*fl[2]+21.78553132241671*fr[1]+21.78553132241671*fl[1]-12.57788237343632*fr[0]+12.57788237343632*fl[0]; 

  incr2[1] = 9.531169953557313*fr[2]+9.531169953557313*fl[2]-2.4609375*fr[1]+2.4609375*fl[1]; 
  incr2[2] = (-36.9140625*fr[2])-36.9140625*fl[2]+9.531169953557313*fr[1]-9.531169953557313*fl[1]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 

} 
