#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 

  incr1[0] = (-0.6708203932499369*fr[7])+0.6708203932499369*fl[7]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[7]-1.161895003862225*fl[7]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = 1.190784930203603*fr[4]+1.190784930203603*fl[4]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = (-2.0625*fr[4])-2.0625*fl[4]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[5] = (-2.0625*fr[5])-2.0625*fl[5]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[6] = 0.9375*fl[6]-0.9375*fr[6]; 
  incr1[7] = (-1.5*fr[7])+1.5*fl[7]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[8] = 0.9375*fl[8]-0.9375*fr[8]; 
  incr1[9] = 0.9375*fl[9]-0.9375*fr[9]; 

  incr2[1] = 0.4236075534914363*fr[7]+0.4236075534914363*fl[7]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.609375*fr[4])+0.609375*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = (-0.609375*fr[5])+0.609375*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-1.640625*fr[7])-1.640625*fl[7]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 

} 
void ConstDiffusionSurf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 

  incr1[0] = (-0.6708203932499369*fr[8])+0.6708203932499369*fl[8]+1.190784930203603*fr[2]+1.190784930203603*fl[2]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.190784930203603*fr[4]+1.190784930203603*fl[4]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.161895003862225*fr[8]-1.161895003862225*fl[8]-2.0625*fr[2]-2.0625*fl[2]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[3] = 1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = (-2.0625*fr[4])-2.0625*fl[4]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[5] = 0.9375*fl[5]-0.9375*fr[5]; 
  incr1[6] = (-2.0625*fr[6])-2.0625*fl[6]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[7] = 0.9375*fl[7]-0.9375*fr[7]; 
  incr1[8] = (-1.5*fr[8])+1.5*fl[8]+2.662676050517599*fr[2]+2.662676050517599*fl[2]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[9] = 0.9375*fl[9]-0.9375*fr[9]; 

  incr2[2] = 0.4236075534914363*fr[8]+0.4236075534914363*fl[8]-0.609375*fr[2]+0.609375*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = (-0.609375*fr[4])+0.609375*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.609375*fr[6])+0.609375*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-1.640625*fr[8])-1.640625*fl[8]+2.360099226595144*fr[2]-2.360099226595144*fl[2]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 

} 
void ConstDiffusionSurf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 

  incr1[0] = (-0.6708203932499369*fr[9])+0.6708203932499369*fl[9]+1.190784930203603*fr[3]+1.190784930203603*fl[3]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.161895003862225*fr[9]-1.161895003862225*fl[9]-2.0625*fr[3]-2.0625*fl[3]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[4] = 0.9375*fl[4]-0.9375*fr[4]; 
  incr1[5] = (-2.0625*fr[5])-2.0625*fl[5]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[6] = (-2.0625*fr[6])-2.0625*fl[6]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[7] = 0.9375*fl[7]-0.9375*fr[7]; 
  incr1[8] = 0.9375*fl[8]-0.9375*fr[8]; 
  incr1[9] = (-1.5*fr[9])+1.5*fl[9]+2.662676050517599*fr[3]+2.662676050517599*fl[3]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 

  incr2[3] = 0.4236075534914363*fr[9]+0.4236075534914363*fl[9]-0.609375*fr[3]+0.609375*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.609375*fr[5])+0.609375*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = (-0.609375*fr[6])+0.609375*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[9] = (-1.640625*fr[9])-1.640625*fl[9]+2.360099226595144*fr[3]-2.360099226595144*fl[3]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  incr1[0] = 6.708203932499369*fr[7]-6.708203932499369*fl[7]-8.11898816047911*fr[1]-8.11898816047911*fl[1]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-11.61895003862225*fr[7])+11.61895003862225*fl[7]+14.0625*fr[1]+14.0625*fl[1]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[2] = (-8.11898816047911*fr[4])-8.11898816047911*fl[4]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-8.11898816047911*fr[5])-8.11898816047911*fl[5]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = 14.0625*fr[4]+14.0625*fl[4]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[5] = 14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[6] = 4.6875*fr[6]-4.6875*fl[6]; 
  incr1[7] = 15.0*fr[7]-15.0*fl[7]-18.15460943534727*fr[1]-18.15460943534727*fl[1]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[8] = 4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = 4.6875*fr[9]-4.6875*fl[9]; 

  incr2[1] = (-2.541645320948617*fr[7])-2.541645320948617*fl[7]+1.40625*fr[1]-1.40625*fl[1]; 
  incr2[4] = 1.40625*fr[4]-1.40625*fl[4]; 
  incr2[5] = 1.40625*fr[5]-1.40625*fl[5]; 
  incr2[7] = 9.84375*fr[7]+9.84375*fl[7]-5.446382830604179*fr[1]+5.446382830604179*fl[1]; 

  incr3[7] = (-4.5*fr[7])+4.5*fl[7]+7.988028151552797*fr[1]+7.988028151552797*fl[1]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += (-1.0*incr2[4]*rdxFnur)-1.0*incr1[4]*rdxFnur; 
  outr[5] += (-1.0*incr2[5]*rdxFnur)-1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += (-1.0*incr3[7]*rdxFnur)-1.0*incr2[7]*rdxFnur-1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul-1.0*incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul-1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul+incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  incr1[0] = 6.708203932499369*fr[8]-6.708203932499369*fl[8]-8.11898816047911*fr[2]-8.11898816047911*fl[2]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-8.11898816047911*fr[4])-8.11898816047911*fl[4]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-11.61895003862225*fr[8])+11.61895003862225*fl[8]+14.0625*fr[2]+14.0625*fl[2]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[3] = (-8.11898816047911*fr[6])-8.11898816047911*fl[6]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = 14.0625*fr[4]+14.0625*fl[4]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[5] = 4.6875*fr[5]-4.6875*fl[5]; 
  incr1[6] = 14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[7] = 4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = 15.0*fr[8]-15.0*fl[8]-18.15460943534727*fr[2]-18.15460943534727*fl[2]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[9] = 4.6875*fr[9]-4.6875*fl[9]; 

  incr2[2] = (-2.541645320948617*fr[8])-2.541645320948617*fl[8]+1.40625*fr[2]-1.40625*fl[2]; 
  incr2[4] = 1.40625*fr[4]-1.40625*fl[4]; 
  incr2[6] = 1.40625*fr[6]-1.40625*fl[6]; 
  incr2[8] = 9.84375*fr[8]+9.84375*fl[8]-5.446382830604179*fr[2]+5.446382830604179*fl[2]; 

  incr3[8] = (-4.5*fr[8])+4.5*fl[8]+7.988028151552797*fr[2]+7.988028151552797*fl[2]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr2[2]*rdxFnur)-1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += (-1.0*incr2[4]*rdxFnur)-1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += (-1.0*incr3[8]*rdxFnur)-1.0*incr2[8]*rdxFnur-1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul-1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr2[8]*rdxFnul+incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  incr1[0] = 6.708203932499369*fr[9]-6.708203932499369*fl[9]-8.11898816047911*fr[3]-8.11898816047911*fl[3]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-8.11898816047911*fr[5])-8.11898816047911*fl[5]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-8.11898816047911*fr[6])-8.11898816047911*fl[6]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-11.61895003862225*fr[9])+11.61895003862225*fl[9]+14.0625*fr[3]+14.0625*fl[3]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[4] = 4.6875*fr[4]-4.6875*fl[4]; 
  incr1[5] = 14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[6] = 14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[7] = 4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = 4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = 15.0*fr[9]-15.0*fl[9]-18.15460943534727*fr[3]-18.15460943534727*fl[3]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 

  incr2[3] = (-2.541645320948617*fr[9])-2.541645320948617*fl[9]+1.40625*fr[3]-1.40625*fl[3]; 
  incr2[5] = 1.40625*fr[5]-1.40625*fl[5]; 
  incr2[6] = 1.40625*fr[6]-1.40625*fl[6]; 
  incr2[9] = 9.84375*fr[9]+9.84375*fl[9]-5.446382830604179*fr[3]+5.446382830604179*fl[3]; 

  incr3[9] = (-4.5*fr[9])+4.5*fl[9]+7.988028151552797*fr[3]+7.988028151552797*fl[3]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr2[3]*rdxFnur)-1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += (-1.0*incr2[5]*rdxFnur)-1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += (-1.0*incr3[9]*rdxFnur)-1.0*incr2[9]*rdxFnur-1.0*incr1[9]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul-1.0*incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr2[9]*rdxFnul+incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  incr1[0] = (-35.21807064562169*fr[7])+35.21807064562169*fl[7]+34.09975027401226*fr[1]+34.09975027401226*fl[1]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 60.9994877027668*fr[7]-60.9994877027668*fl[7]-59.0625*fr[1]-59.0625*fl[1]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[2] = 34.09975027401226*fr[4]+34.09975027401226*fl[4]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = (-59.0625*fr[4])-59.0625*fl[4]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[5] = (-59.0625*fr[5])-59.0625*fl[5]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[6] = 19.6875*fl[6]-19.6875*fr[6]; 
  incr1[7] = (-78.75*fr[7])+78.75*fl[7]+76.2493596284585*fr[1]+76.2493596284585*fl[1]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[8] = 19.6875*fl[8]-19.6875*fr[8]; 
  incr1[9] = 19.6875*fl[9]-19.6875*fr[9]; 

  incr2[1] = (-9.531169953557313*fr[7])-9.531169953557313*fl[7]+2.4609375*fr[1]-2.4609375*fl[1]; 
  incr2[4] = 2.4609375*fr[4]-2.4609375*fl[4]; 
  incr2[5] = 2.4609375*fr[5]-2.4609375*fl[5]; 
  incr2[7] = 36.9140625*fr[7]+36.9140625*fl[7]-9.531169953557313*fr[1]+9.531169953557313*fl[1]; 

  incr3[7] = (-45.0*fr[7])+45.0*fl[7]+54.4638283060418*fr[1]+54.4638283060418*fl[1]-31.4447059335908*fr[0]+31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr3[7]*rdxFnur+incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul+incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul+incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul+incr1[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  incr1[0] = (-35.21807064562169*fr[8])+35.21807064562169*fl[8]+34.09975027401226*fr[2]+34.09975027401226*fl[2]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 34.09975027401226*fr[4]+34.09975027401226*fl[4]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 60.9994877027668*fr[8]-60.9994877027668*fl[8]-59.0625*fr[2]-59.0625*fl[2]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[3] = 34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = (-59.0625*fr[4])-59.0625*fl[4]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[5] = 19.6875*fl[5]-19.6875*fr[5]; 
  incr1[6] = (-59.0625*fr[6])-59.0625*fl[6]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[7] = 19.6875*fl[7]-19.6875*fr[7]; 
  incr1[8] = (-78.75*fr[8])+78.75*fl[8]+76.2493596284585*fr[2]+76.2493596284585*fl[2]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[9] = 19.6875*fl[9]-19.6875*fr[9]; 

  incr2[2] = (-9.531169953557313*fr[8])-9.531169953557313*fl[8]+2.4609375*fr[2]-2.4609375*fl[2]; 
  incr2[4] = 2.4609375*fr[4]-2.4609375*fl[4]; 
  incr2[6] = 2.4609375*fr[6]-2.4609375*fl[6]; 
  incr2[8] = 36.9140625*fr[8]+36.9140625*fl[8]-9.531169953557313*fr[2]+9.531169953557313*fl[2]; 

  incr3[8] = (-45.0*fr[8])+45.0*fl[8]+54.4638283060418*fr[2]+54.4638283060418*fl[2]-31.4447059335908*fr[0]+31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr3[8]*rdxFnur+incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul+incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul+incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul+incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  incr1[0] = (-35.21807064562169*fr[9])+35.21807064562169*fl[9]+34.09975027401226*fr[3]+34.09975027401226*fl[3]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 60.9994877027668*fr[9]-60.9994877027668*fl[9]-59.0625*fr[3]-59.0625*fl[3]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[4] = 19.6875*fl[4]-19.6875*fr[4]; 
  incr1[5] = (-59.0625*fr[5])-59.0625*fl[5]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[6] = (-59.0625*fr[6])-59.0625*fl[6]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[7] = 19.6875*fl[7]-19.6875*fr[7]; 
  incr1[8] = 19.6875*fl[8]-19.6875*fr[8]; 
  incr1[9] = (-78.75*fr[9])+78.75*fl[9]+76.2493596284585*fr[3]+76.2493596284585*fl[3]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 

  incr2[3] = (-9.531169953557313*fr[9])-9.531169953557313*fl[9]+2.4609375*fr[3]-2.4609375*fl[3]; 
  incr2[5] = 2.4609375*fr[5]-2.4609375*fl[5]; 
  incr2[6] = 2.4609375*fr[6]-2.4609375*fl[6]; 
  incr2[9] = 36.9140625*fr[9]+36.9140625*fl[9]-9.531169953557313*fr[3]+9.531169953557313*fl[3]; 

  incr3[9] = (-45.0*fr[9])+45.0*fl[9]+54.4638283060418*fr[3]+54.4638283060418*fl[3]-31.4447059335908*fr[0]+31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr3[9]*rdxFnur+incr2[9]*rdxFnur+incr1[9]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul+incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul+incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul+incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 

} 
