#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 

  incr1[0] = (-0.6708203932499369*fr[11])+0.6708203932499369*fl[11]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[11]-1.161895003862225*fl[11]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = 1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = 1.190784930203603*fr[8]+1.190784930203603*fl[8]-0.9375*fr[4]+0.9375*fl[4]; 
  incr1[5] = (-2.0625*fr[5])-2.0625*fl[5]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[6] = (-2.0625*fr[6])-2.0625*fl[6]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[7] = 0.9375*fl[7]-0.9375*fr[7]; 
  incr1[8] = (-2.0625*fr[8])-2.0625*fl[8]+1.623797632095822*fr[4]-1.623797632095822*fl[4]; 
  incr1[9] = 0.9375*fl[9]-0.9375*fr[9]; 
  incr1[10] = 0.9375*fl[10]-0.9375*fr[10]; 
  incr1[11] = (-1.5*fr[11])+1.5*fl[11]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[12] = 0.9375*fl[12]-0.9375*fr[12]; 
  incr1[13] = 0.9375*fl[13]-0.9375*fr[13]; 
  incr1[14] = 0.9375*fl[14]-0.9375*fr[14]; 

  incr2[1] = 0.4236075534914363*fr[11]+0.4236075534914363*fl[11]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.609375*fr[5])+0.609375*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.609375*fr[6])+0.609375*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.609375*fr[8])+0.609375*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-1.640625*fr[11])-1.640625*fl[11]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstDiffusionSurf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 

  incr1[0] = (-0.6708203932499369*fr[12])+0.6708203932499369*fl[12]+1.190784930203603*fr[2]+1.190784930203603*fl[2]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.161895003862225*fr[12]-1.161895003862225*fl[12]-2.0625*fr[2]-2.0625*fl[2]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[3] = 1.190784930203603*fr[7]+1.190784930203603*fl[7]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = 1.190784930203603*fr[9]+1.190784930203603*fl[9]-0.9375*fr[4]+0.9375*fl[4]; 
  incr1[5] = (-2.0625*fr[5])-2.0625*fl[5]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[6] = 0.9375*fl[6]-0.9375*fr[6]; 
  incr1[7] = (-2.0625*fr[7])-2.0625*fl[7]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[8] = 0.9375*fl[8]-0.9375*fr[8]; 
  incr1[9] = (-2.0625*fr[9])-2.0625*fl[9]+1.623797632095822*fr[4]-1.623797632095822*fl[4]; 
  incr1[10] = 0.9375*fl[10]-0.9375*fr[10]; 
  incr1[11] = 0.9375*fl[11]-0.9375*fr[11]; 
  incr1[12] = (-1.5*fr[12])+1.5*fl[12]+2.662676050517599*fr[2]+2.662676050517599*fl[2]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[13] = 0.9375*fl[13]-0.9375*fr[13]; 
  incr1[14] = 0.9375*fl[14]-0.9375*fr[14]; 

  incr2[2] = 0.4236075534914363*fr[12]+0.4236075534914363*fl[12]-0.609375*fr[2]+0.609375*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.609375*fr[5])+0.609375*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.609375*fr[7])+0.609375*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.609375*fr[9])+0.609375*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[12] = (-1.640625*fr[12])-1.640625*fl[12]+2.360099226595144*fr[2]-2.360099226595144*fl[2]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstDiffusionSurf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 

  incr1[0] = (-0.6708203932499369*fr[13])+0.6708203932499369*fl[13]+1.190784930203603*fr[3]+1.190784930203603*fl[3]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.190784930203603*fr[7]+1.190784930203603*fl[7]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.161895003862225*fr[13]-1.161895003862225*fl[13]-2.0625*fr[3]-2.0625*fl[3]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[4] = 1.190784930203603*fr[10]+1.190784930203603*fl[10]-0.9375*fr[4]+0.9375*fl[4]; 
  incr1[5] = 0.9375*fl[5]-0.9375*fr[5]; 
  incr1[6] = (-2.0625*fr[6])-2.0625*fl[6]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[7] = (-2.0625*fr[7])-2.0625*fl[7]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[8] = 0.9375*fl[8]-0.9375*fr[8]; 
  incr1[9] = 0.9375*fl[9]-0.9375*fr[9]; 
  incr1[10] = (-2.0625*fr[10])-2.0625*fl[10]+1.623797632095822*fr[4]-1.623797632095822*fl[4]; 
  incr1[11] = 0.9375*fl[11]-0.9375*fr[11]; 
  incr1[12] = 0.9375*fl[12]-0.9375*fr[12]; 
  incr1[13] = (-1.5*fr[13])+1.5*fl[13]+2.662676050517599*fr[3]+2.662676050517599*fl[3]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[14] = 0.9375*fl[14]-0.9375*fr[14]; 

  incr2[3] = 0.4236075534914363*fr[13]+0.4236075534914363*fl[13]-0.609375*fr[3]+0.609375*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.609375*fr[6])+0.609375*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.609375*fr[7])+0.609375*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.609375*fr[10])+0.609375*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[13] = (-1.640625*fr[13])-1.640625*fl[13]+2.360099226595144*fr[3]-2.360099226595144*fl[3]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstDiffusionSurf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 

  incr1[0] = (-0.6708203932499369*fr[14])+0.6708203932499369*fl[14]+1.190784930203603*fr[4]+1.190784930203603*fl[4]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.190784930203603*fr[8]+1.190784930203603*fl[8]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.190784930203603*fr[9]+1.190784930203603*fl[9]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.190784930203603*fr[10]+1.190784930203603*fl[10]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = 1.161895003862225*fr[14]-1.161895003862225*fl[14]-2.0625*fr[4]-2.0625*fl[4]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[5] = 0.9375*fl[5]-0.9375*fr[5]; 
  incr1[6] = 0.9375*fl[6]-0.9375*fr[6]; 
  incr1[7] = 0.9375*fl[7]-0.9375*fr[7]; 
  incr1[8] = (-2.0625*fr[8])-2.0625*fl[8]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[9] = (-2.0625*fr[9])-2.0625*fl[9]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[10] = (-2.0625*fr[10])-2.0625*fl[10]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[11] = 0.9375*fl[11]-0.9375*fr[11]; 
  incr1[12] = 0.9375*fl[12]-0.9375*fr[12]; 
  incr1[13] = 0.9375*fl[13]-0.9375*fr[13]; 
  incr1[14] = (-1.5*fr[14])+1.5*fl[14]+2.662676050517599*fr[4]+2.662676050517599*fl[4]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 

  incr2[4] = 0.4236075534914363*fr[14]+0.4236075534914363*fl[14]-0.609375*fr[4]+0.609375*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.609375*fr[8])+0.609375*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.609375*fr[9])+0.609375*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.609375*fr[10])+0.609375*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[14] = (-1.640625*fr[14])-1.640625*fl[14]+2.360099226595144*fr[4]-2.360099226595144*fl[4]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  incr1[0] = 6.708203932499369*fr[11]-6.708203932499369*fl[11]-8.11898816047911*fr[1]-8.11898816047911*fl[1]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-11.61895003862225*fr[11])+11.61895003862225*fl[11]+14.0625*fr[1]+14.0625*fl[1]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[2] = (-8.11898816047911*fr[5])-8.11898816047911*fl[5]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-8.11898816047911*fr[6])-8.11898816047911*fl[6]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = (-8.11898816047911*fr[8])-8.11898816047911*fl[8]+4.6875*fr[4]-4.6875*fl[4]; 
  incr1[5] = 14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[6] = 14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[7] = 4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = 14.0625*fr[8]+14.0625*fl[8]-8.11898816047911*fr[4]+8.11898816047911*fl[4]; 
  incr1[9] = 4.6875*fr[9]-4.6875*fl[9]; 
  incr1[10] = 4.6875*fr[10]-4.6875*fl[10]; 
  incr1[11] = 15.0*fr[11]-15.0*fl[11]-18.15460943534727*fr[1]-18.15460943534727*fl[1]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[12] = 4.6875*fr[12]-4.6875*fl[12]; 
  incr1[13] = 4.6875*fr[13]-4.6875*fl[13]; 
  incr1[14] = 4.6875*fr[14]-4.6875*fl[14]; 

  incr2[1] = (-2.541645320948617*fr[11])-2.541645320948617*fl[11]+1.40625*fr[1]-1.40625*fl[1]; 
  incr2[5] = 1.40625*fr[5]-1.40625*fl[5]; 
  incr2[6] = 1.40625*fr[6]-1.40625*fl[6]; 
  incr2[8] = 1.40625*fr[8]-1.40625*fl[8]; 
  incr2[11] = 9.84375*fr[11]+9.84375*fl[11]-5.446382830604179*fr[1]+5.446382830604179*fl[1]; 

  incr3[11] = (-4.5*fr[11])+4.5*fl[11]+7.988028151552797*fr[1]+7.988028151552797*fl[1]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += (-1.0*incr2[5]*rdxFnur)-1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += (-1.0*incr2[8]*rdxFnur)-1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 
  outr[10] += -1.0*incr1[10]*rdxFnur; 
  outr[11] += (-1.0*incr3[11]*rdxFnur)-1.0*incr2[11]*rdxFnur-1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul-1.0*incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr2[11]*rdxFnul+incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  incr1[0] = 6.708203932499369*fr[12]-6.708203932499369*fl[12]-8.11898816047911*fr[2]-8.11898816047911*fl[2]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-8.11898816047911*fr[5])-8.11898816047911*fl[5]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-11.61895003862225*fr[12])+11.61895003862225*fl[12]+14.0625*fr[2]+14.0625*fl[2]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[3] = (-8.11898816047911*fr[7])-8.11898816047911*fl[7]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = (-8.11898816047911*fr[9])-8.11898816047911*fl[9]+4.6875*fr[4]-4.6875*fl[4]; 
  incr1[5] = 14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[6] = 4.6875*fr[6]-4.6875*fl[6]; 
  incr1[7] = 14.0625*fr[7]+14.0625*fl[7]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[8] = 4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = 14.0625*fr[9]+14.0625*fl[9]-8.11898816047911*fr[4]+8.11898816047911*fl[4]; 
  incr1[10] = 4.6875*fr[10]-4.6875*fl[10]; 
  incr1[11] = 4.6875*fr[11]-4.6875*fl[11]; 
  incr1[12] = 15.0*fr[12]-15.0*fl[12]-18.15460943534727*fr[2]-18.15460943534727*fl[2]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[13] = 4.6875*fr[13]-4.6875*fl[13]; 
  incr1[14] = 4.6875*fr[14]-4.6875*fl[14]; 

  incr2[2] = (-2.541645320948617*fr[12])-2.541645320948617*fl[12]+1.40625*fr[2]-1.40625*fl[2]; 
  incr2[5] = 1.40625*fr[5]-1.40625*fl[5]; 
  incr2[7] = 1.40625*fr[7]-1.40625*fl[7]; 
  incr2[9] = 1.40625*fr[9]-1.40625*fl[9]; 
  incr2[12] = 9.84375*fr[12]+9.84375*fl[12]-5.446382830604179*fr[2]+5.446382830604179*fl[2]; 

  incr3[12] = (-4.5*fr[12])+4.5*fl[12]+7.988028151552797*fr[2]+7.988028151552797*fl[2]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr2[2]*rdxFnur)-1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += (-1.0*incr2[5]*rdxFnur)-1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += (-1.0*incr2[7]*rdxFnur)-1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += (-1.0*incr2[9]*rdxFnur)-1.0*incr1[9]*rdxFnur; 
  outr[10] += -1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr3[12]*rdxFnur)-1.0*incr2[12]*rdxFnur-1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul-1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr2[12]*rdxFnul+incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  incr1[0] = 6.708203932499369*fr[13]-6.708203932499369*fl[13]-8.11898816047911*fr[3]-8.11898816047911*fl[3]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-8.11898816047911*fr[6])-8.11898816047911*fl[6]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-8.11898816047911*fr[7])-8.11898816047911*fl[7]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-11.61895003862225*fr[13])+11.61895003862225*fl[13]+14.0625*fr[3]+14.0625*fl[3]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[4] = (-8.11898816047911*fr[10])-8.11898816047911*fl[10]+4.6875*fr[4]-4.6875*fl[4]; 
  incr1[5] = 4.6875*fr[5]-4.6875*fl[5]; 
  incr1[6] = 14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[7] = 14.0625*fr[7]+14.0625*fl[7]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[8] = 4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = 4.6875*fr[9]-4.6875*fl[9]; 
  incr1[10] = 14.0625*fr[10]+14.0625*fl[10]-8.11898816047911*fr[4]+8.11898816047911*fl[4]; 
  incr1[11] = 4.6875*fr[11]-4.6875*fl[11]; 
  incr1[12] = 4.6875*fr[12]-4.6875*fl[12]; 
  incr1[13] = 15.0*fr[13]-15.0*fl[13]-18.15460943534727*fr[3]-18.15460943534727*fl[3]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[14] = 4.6875*fr[14]-4.6875*fl[14]; 

  incr2[3] = (-2.541645320948617*fr[13])-2.541645320948617*fl[13]+1.40625*fr[3]-1.40625*fl[3]; 
  incr2[6] = 1.40625*fr[6]-1.40625*fl[6]; 
  incr2[7] = 1.40625*fr[7]-1.40625*fl[7]; 
  incr2[10] = 1.40625*fr[10]-1.40625*fl[10]; 
  incr2[13] = 9.84375*fr[13]+9.84375*fl[13]-5.446382830604179*fr[3]+5.446382830604179*fl[3]; 

  incr3[13] = (-4.5*fr[13])+4.5*fl[13]+7.988028151552797*fr[3]+7.988028151552797*fl[3]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr2[3]*rdxFnur)-1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += (-1.0*incr2[7]*rdxFnur)-1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur-1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul+incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  incr1[0] = 6.708203932499369*fr[14]-6.708203932499369*fl[14]-8.11898816047911*fr[4]-8.11898816047911*fl[4]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-8.11898816047911*fr[8])-8.11898816047911*fl[8]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-8.11898816047911*fr[9])-8.11898816047911*fl[9]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-8.11898816047911*fr[10])-8.11898816047911*fl[10]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = (-11.61895003862225*fr[14])+11.61895003862225*fl[14]+14.0625*fr[4]+14.0625*fl[4]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[5] = 4.6875*fr[5]-4.6875*fl[5]; 
  incr1[6] = 4.6875*fr[6]-4.6875*fl[6]; 
  incr1[7] = 4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = 14.0625*fr[8]+14.0625*fl[8]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[9] = 14.0625*fr[9]+14.0625*fl[9]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[10] = 14.0625*fr[10]+14.0625*fl[10]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[11] = 4.6875*fr[11]-4.6875*fl[11]; 
  incr1[12] = 4.6875*fr[12]-4.6875*fl[12]; 
  incr1[13] = 4.6875*fr[13]-4.6875*fl[13]; 
  incr1[14] = 15.0*fr[14]-15.0*fl[14]-18.15460943534727*fr[4]-18.15460943534727*fl[4]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 

  incr2[4] = (-2.541645320948617*fr[14])-2.541645320948617*fl[14]+1.40625*fr[4]-1.40625*fl[4]; 
  incr2[8] = 1.40625*fr[8]-1.40625*fl[8]; 
  incr2[9] = 1.40625*fr[9]-1.40625*fl[9]; 
  incr2[10] = 1.40625*fr[10]-1.40625*fl[10]; 
  incr2[14] = 9.84375*fr[14]+9.84375*fl[14]-5.446382830604179*fr[4]+5.446382830604179*fl[4]; 

  incr3[14] = (-4.5*fr[14])+4.5*fl[14]+7.988028151552797*fr[4]+7.988028151552797*fl[4]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 


  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += (-1.0*incr2[4]*rdxFnur)-1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += (-1.0*incr2[8]*rdxFnur)-1.0*incr1[8]*rdxFnur; 
  outr[9] += (-1.0*incr2[9]*rdxFnur)-1.0*incr1[9]*rdxFnur; 
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur-1.0*incr1[14]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul-1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul+incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  incr1[0] = (-35.21807064562169*fr[11])+35.21807064562169*fl[11]+34.09975027401226*fr[1]+34.09975027401226*fl[1]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 60.9994877027668*fr[11]-60.9994877027668*fl[11]-59.0625*fr[1]-59.0625*fl[1]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[2] = 34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = 34.09975027401226*fr[8]+34.09975027401226*fl[8]-19.6875*fr[4]+19.6875*fl[4]; 
  incr1[5] = (-59.0625*fr[5])-59.0625*fl[5]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[6] = (-59.0625*fr[6])-59.0625*fl[6]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[7] = 19.6875*fl[7]-19.6875*fr[7]; 
  incr1[8] = (-59.0625*fr[8])-59.0625*fl[8]+34.09975027401226*fr[4]-34.09975027401226*fl[4]; 
  incr1[9] = 19.6875*fl[9]-19.6875*fr[9]; 
  incr1[10] = 19.6875*fl[10]-19.6875*fr[10]; 
  incr1[11] = (-78.75*fr[11])+78.75*fl[11]+76.2493596284585*fr[1]+76.2493596284585*fl[1]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[12] = 19.6875*fl[12]-19.6875*fr[12]; 
  incr1[13] = 19.6875*fl[13]-19.6875*fr[13]; 
  incr1[14] = 19.6875*fl[14]-19.6875*fr[14]; 

  incr2[1] = 9.531169953557313*fr[11]+9.531169953557313*fl[11]-2.4609375*fr[1]+2.4609375*fl[1]; 
  incr2[5] = 2.4609375*fl[5]-2.4609375*fr[5]; 
  incr2[6] = 2.4609375*fl[6]-2.4609375*fr[6]; 
  incr2[8] = 2.4609375*fl[8]-2.4609375*fr[8]; 
  incr2[11] = (-36.9140625*fr[11])-36.9140625*fl[11]+9.531169953557313*fr[1]-9.531169953557313*fl[1]; 

  incr3[11] = 45.0*fr[11]-45.0*fl[11]-54.4638283060418*fr[1]-54.4638283060418*fl[1]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr3[11]*rdxFnur+incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += (-1.0*incr3[11]*rdxFnul)+incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  incr1[0] = (-35.21807064562169*fr[12])+35.21807064562169*fl[12]+34.09975027401226*fr[2]+34.09975027401226*fl[2]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 60.9994877027668*fr[12]-60.9994877027668*fl[12]-59.0625*fr[2]-59.0625*fl[2]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[3] = 34.09975027401226*fr[7]+34.09975027401226*fl[7]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = 34.09975027401226*fr[9]+34.09975027401226*fl[9]-19.6875*fr[4]+19.6875*fl[4]; 
  incr1[5] = (-59.0625*fr[5])-59.0625*fl[5]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[6] = 19.6875*fl[6]-19.6875*fr[6]; 
  incr1[7] = (-59.0625*fr[7])-59.0625*fl[7]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[8] = 19.6875*fl[8]-19.6875*fr[8]; 
  incr1[9] = (-59.0625*fr[9])-59.0625*fl[9]+34.09975027401226*fr[4]-34.09975027401226*fl[4]; 
  incr1[10] = 19.6875*fl[10]-19.6875*fr[10]; 
  incr1[11] = 19.6875*fl[11]-19.6875*fr[11]; 
  incr1[12] = (-78.75*fr[12])+78.75*fl[12]+76.2493596284585*fr[2]+76.2493596284585*fl[2]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[13] = 19.6875*fl[13]-19.6875*fr[13]; 
  incr1[14] = 19.6875*fl[14]-19.6875*fr[14]; 

  incr2[2] = 9.531169953557313*fr[12]+9.531169953557313*fl[12]-2.4609375*fr[2]+2.4609375*fl[2]; 
  incr2[5] = 2.4609375*fl[5]-2.4609375*fr[5]; 
  incr2[7] = 2.4609375*fl[7]-2.4609375*fr[7]; 
  incr2[9] = 2.4609375*fl[9]-2.4609375*fr[9]; 
  incr2[12] = (-36.9140625*fr[12])-36.9140625*fl[12]+9.531169953557313*fr[2]-9.531169953557313*fl[2]; 

  incr3[12] = 45.0*fr[12]-45.0*fl[12]-54.4638283060418*fr[2]-54.4638283060418*fl[2]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr3[12]*rdxFnur+incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += (-1.0*incr3[12]*rdxFnul)+incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  incr1[0] = (-35.21807064562169*fr[13])+35.21807064562169*fl[13]+34.09975027401226*fr[3]+34.09975027401226*fl[3]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 34.09975027401226*fr[7]+34.09975027401226*fl[7]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 60.9994877027668*fr[13]-60.9994877027668*fl[13]-59.0625*fr[3]-59.0625*fl[3]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[4] = 34.09975027401226*fr[10]+34.09975027401226*fl[10]-19.6875*fr[4]+19.6875*fl[4]; 
  incr1[5] = 19.6875*fl[5]-19.6875*fr[5]; 
  incr1[6] = (-59.0625*fr[6])-59.0625*fl[6]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[7] = (-59.0625*fr[7])-59.0625*fl[7]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[8] = 19.6875*fl[8]-19.6875*fr[8]; 
  incr1[9] = 19.6875*fl[9]-19.6875*fr[9]; 
  incr1[10] = (-59.0625*fr[10])-59.0625*fl[10]+34.09975027401226*fr[4]-34.09975027401226*fl[4]; 
  incr1[11] = 19.6875*fl[11]-19.6875*fr[11]; 
  incr1[12] = 19.6875*fl[12]-19.6875*fr[12]; 
  incr1[13] = (-78.75*fr[13])+78.75*fl[13]+76.2493596284585*fr[3]+76.2493596284585*fl[3]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[14] = 19.6875*fl[14]-19.6875*fr[14]; 

  incr2[3] = 9.531169953557313*fr[13]+9.531169953557313*fl[13]-2.4609375*fr[3]+2.4609375*fl[3]; 
  incr2[6] = 2.4609375*fl[6]-2.4609375*fr[6]; 
  incr2[7] = 2.4609375*fl[7]-2.4609375*fr[7]; 
  incr2[10] = 2.4609375*fl[10]-2.4609375*fr[10]; 
  incr2[13] = (-36.9140625*fr[13])-36.9140625*fl[13]+9.531169953557313*fr[3]-9.531169953557313*fl[3]; 

  incr3[13] = 45.0*fr[13]-45.0*fl[13]-54.4638283060418*fr[3]-54.4638283060418*fl[3]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr3[13]*rdxFnur+incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += (-1.0*incr3[13]*rdxFnul)+incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 64.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  incr1[0] = (-35.21807064562169*fr[14])+35.21807064562169*fl[14]+34.09975027401226*fr[4]+34.09975027401226*fl[4]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 34.09975027401226*fr[8]+34.09975027401226*fl[8]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 34.09975027401226*fr[9]+34.09975027401226*fl[9]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 34.09975027401226*fr[10]+34.09975027401226*fl[10]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = 60.9994877027668*fr[14]-60.9994877027668*fl[14]-59.0625*fr[4]-59.0625*fl[4]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[5] = 19.6875*fl[5]-19.6875*fr[5]; 
  incr1[6] = 19.6875*fl[6]-19.6875*fr[6]; 
  incr1[7] = 19.6875*fl[7]-19.6875*fr[7]; 
  incr1[8] = (-59.0625*fr[8])-59.0625*fl[8]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[9] = (-59.0625*fr[9])-59.0625*fl[9]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[10] = (-59.0625*fr[10])-59.0625*fl[10]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[11] = 19.6875*fl[11]-19.6875*fr[11]; 
  incr1[12] = 19.6875*fl[12]-19.6875*fr[12]; 
  incr1[13] = 19.6875*fl[13]-19.6875*fr[13]; 
  incr1[14] = (-78.75*fr[14])+78.75*fl[14]+76.2493596284585*fr[4]+76.2493596284585*fl[4]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 

  incr2[4] = 9.531169953557313*fr[14]+9.531169953557313*fl[14]-2.4609375*fr[4]+2.4609375*fl[4]; 
  incr2[8] = 2.4609375*fl[8]-2.4609375*fr[8]; 
  incr2[9] = 2.4609375*fl[9]-2.4609375*fr[9]; 
  incr2[10] = 2.4609375*fl[10]-2.4609375*fr[10]; 
  incr2[14] = (-36.9140625*fr[14])-36.9140625*fl[14]+9.531169953557313*fr[4]-9.531169953557313*fl[4]; 

  incr3[14] = 45.0*fr[14]-45.0*fl[14]-54.4638283060418*fr[4]-54.4638283060418*fl[4]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 




  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr3[14]*rdxFnur+incr2[14]*rdxFnur+incr1[14]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += (-1.0*incr3[14]*rdxFnul)+incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 

} 
