#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf3xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxvr[0]*dxvr[0]); 

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
void ConstDiffusionSurf3xSer_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[20]; 
  double incr2[20]; 

  incr1[0] = (-0.6708203932499369*fr[7])+0.6708203932499369*fl[7]+1.190784930203603*fr[1]+1.190784930203603*fl[1]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = 1.161895003862225*fr[7]-1.161895003862225*fl[7]-2.0625*fr[1]-2.0625*fl[1]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[2] = (-0.6708203932499369*fr[11])+0.6708203932499369*fl[11]+1.190784930203603*fr[4]+1.190784930203603*fl[4]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = (-0.6708203932499369*fr[13])+0.6708203932499369*fl[13]+1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = 1.161895003862225*fr[11]-1.161895003862225*fl[11]-2.0625*fr[4]-2.0625*fl[4]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[5] = 1.161895003862225*fr[13]-1.161895003862225*fl[13]-2.0625*fr[5]-2.0625*fl[5]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[6] = (-0.6708203932499369*fr[17])+0.6708203932499369*fl[17]+1.190784930203603*fr[10]+1.190784930203603*fl[10]-0.9375*fr[6]+0.9375*fl[6]; 
  incr1[7] = (-1.5*fr[7])+1.5*fl[7]+2.662676050517599*fr[1]+2.662676050517599*fl[1]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[8] = 1.190784930203603*fr[12]+1.190784930203603*fl[12]-0.9375*fr[8]+0.9375*fl[8]; 
  incr1[9] = 1.190784930203603*fr[15]+1.190784930203603*fl[15]-0.9375*fr[9]+0.9375*fl[9]; 
  incr1[10] = 1.161895003862225*fr[17]-1.161895003862225*fl[17]-2.0625*fr[10]-2.0625*fl[10]+1.623797632095822*fr[6]-1.623797632095822*fl[6]; 
  incr1[11] = (-1.5*fr[11])+1.5*fl[11]+2.662676050517599*fr[4]+2.662676050517599*fl[4]-2.096313728906053*fr[2]+2.096313728906053*fl[2]; 
  incr1[12] = (-2.0625*fr[12])-2.0625*fl[12]+1.623797632095823*fr[8]-1.623797632095823*fl[8]; 
  incr1[13] = (-1.5*fr[13])+1.5*fl[13]+2.662676050517599*fr[5]+2.662676050517599*fl[5]-2.096313728906053*fr[3]+2.096313728906053*fl[3]; 
  incr1[14] = 1.190784930203603*fr[18]+1.190784930203603*fl[18]-0.9375*fr[14]+0.9375*fl[14]; 
  incr1[15] = (-2.0625*fr[15])-2.0625*fl[15]+1.623797632095823*fr[9]-1.623797632095823*fl[9]; 
  incr1[16] = 1.190784930203603*fr[19]+1.190784930203603*fl[19]-0.9375*fr[16]+0.9375*fl[16]; 
  incr1[17] = (-1.5*fr[17])+1.5*fl[17]+2.662676050517599*fr[10]+2.662676050517599*fl[10]-2.096313728906053*fr[6]+2.096313728906053*fl[6]; 
  incr1[18] = (-2.0625*fr[18])-2.0625*fl[18]+1.623797632095823*fr[14]-1.623797632095823*fl[14]; 
  incr1[19] = (-2.0625*fr[19])-2.0625*fl[19]+1.623797632095823*fr[16]-1.623797632095823*fl[16]; 

  incr2[1] = 0.4236075534914363*fr[7]+0.4236075534914363*fl[7]-0.609375*fr[1]+0.609375*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = 0.4236075534914363*fr[11]+0.4236075534914363*fl[11]-0.609375*fr[4]+0.609375*fl[4]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[5] = 0.4236075534914363*fr[13]+0.4236075534914363*fl[13]-0.609375*fr[5]+0.609375*fl[5]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[7] = (-1.640625*fr[7])-1.640625*fl[7]+2.360099226595144*fr[1]-2.360099226595144*fl[1]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 
  incr2[10] = 0.4236075534914363*fr[17]+0.4236075534914363*fl[17]-0.609375*fr[10]+0.609375*fl[10]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[11] = (-1.640625*fr[11])-1.640625*fl[11]+2.360099226595145*fr[4]-2.360099226595145*fl[4]-1.677050983124842*fr[2]-1.677050983124842*fl[2]; 
  incr2[12] = (-0.609375*fr[12])+0.609375*fl[12]+0.4330127018922194*fr[8]+0.4330127018922194*fl[8]; 
  incr2[13] = (-1.640625*fr[13])-1.640625*fl[13]+2.360099226595145*fr[5]-2.360099226595145*fl[5]-1.677050983124842*fr[3]-1.677050983124842*fl[3]; 
  incr2[15] = (-0.609375*fr[15])+0.609375*fl[15]+0.4330127018922194*fr[9]+0.4330127018922194*fl[9]; 
  incr2[17] = (-1.640625*fr[17])-1.640625*fl[17]+2.360099226595144*fr[10]-2.360099226595144*fl[10]-1.677050983124842*fr[6]-1.677050983124842*fl[6]; 
  incr2[18] = (-0.609375*fr[18])+0.609375*fl[18]+0.4330127018922194*fr[14]+0.4330127018922194*fl[14]; 
  incr2[19] = (-0.609375*fr[19])+0.609375*fl[19]+0.4330127018922194*fr[16]+0.4330127018922194*fl[16]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 
  outr[16] += incr1[16]*rdxSq4nur; 
  outr[17] += incr2[17]*rdxSq4nur+incr1[17]*rdxSq4nur; 
  outr[18] += incr2[18]*rdxSq4nur+incr1[18]*rdxSq4nur; 
  outr[19] += incr2[19]*rdxSq4nur+incr1[19]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr2[7]*rdxSq4nul-1.0*incr1[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr2[11]*rdxSq4nul-1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr2[13]*rdxSq4nul-1.0*incr1[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 
  outl[16] += -1.0*incr1[16]*rdxSq4nul; 
  outl[17] += incr2[17]*rdxSq4nul-1.0*incr1[17]*rdxSq4nul; 
  outl[18] += incr1[18]*rdxSq4nul-1.0*incr2[18]*rdxSq4nul; 
  outl[19] += incr1[19]*rdxSq4nul-1.0*incr2[19]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxvr[1]*dxvr[1]); 

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
void ConstDiffusionSurf3xSer_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[20]; 
  double incr2[20]; 

  incr1[0] = (-0.6708203932499369*fr[8])+0.6708203932499369*fl[8]+1.190784930203603*fr[2]+1.190784930203603*fl[2]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = (-0.6708203932499369*fr[12])+0.6708203932499369*fl[12]+1.190784930203603*fr[4]+1.190784930203603*fl[4]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = 1.161895003862225*fr[8]-1.161895003862225*fl[8]-2.0625*fr[2]-2.0625*fl[2]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[3] = (-0.6708203932499369*fr[14])+0.6708203932499369*fl[14]+1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[3]+0.9375*fl[3]; 
  incr1[4] = 1.161895003862225*fr[12]-1.161895003862225*fl[12]-2.0625*fr[4]-2.0625*fl[4]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[5] = (-0.6708203932499369*fr[18])+0.6708203932499369*fl[18]+1.190784930203603*fr[10]+1.190784930203603*fl[10]-0.9375*fr[5]+0.9375*fl[5]; 
  incr1[6] = 1.161895003862225*fr[14]-1.161895003862225*fl[14]-2.0625*fr[6]-2.0625*fl[6]+1.623797632095822*fr[3]-1.623797632095822*fl[3]; 
  incr1[7] = 1.190784930203603*fr[11]+1.190784930203603*fl[11]-0.9375*fr[7]+0.9375*fl[7]; 
  incr1[8] = (-1.5*fr[8])+1.5*fl[8]+2.662676050517599*fr[2]+2.662676050517599*fl[2]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[9] = 1.190784930203603*fr[16]+1.190784930203603*fl[16]-0.9375*fr[9]+0.9375*fl[9]; 
  incr1[10] = 1.161895003862225*fr[18]-1.161895003862225*fl[18]-2.0625*fr[10]-2.0625*fl[10]+1.623797632095822*fr[5]-1.623797632095822*fl[5]; 
  incr1[11] = (-2.0625*fr[11])-2.0625*fl[11]+1.623797632095823*fr[7]-1.623797632095823*fl[7]; 
  incr1[12] = (-1.5*fr[12])+1.5*fl[12]+2.662676050517599*fr[4]+2.662676050517599*fl[4]-2.096313728906053*fr[1]+2.096313728906053*fl[1]; 
  incr1[13] = 1.190784930203603*fr[17]+1.190784930203603*fl[17]-0.9375*fr[13]+0.9375*fl[13]; 
  incr1[14] = (-1.5*fr[14])+1.5*fl[14]+2.662676050517599*fr[6]+2.662676050517599*fl[6]-2.096313728906053*fr[3]+2.096313728906053*fl[3]; 
  incr1[15] = 1.190784930203603*fr[19]+1.190784930203603*fl[19]-0.9375*fr[15]+0.9375*fl[15]; 
  incr1[16] = (-2.0625*fr[16])-2.0625*fl[16]+1.623797632095823*fr[9]-1.623797632095823*fl[9]; 
  incr1[17] = (-2.0625*fr[17])-2.0625*fl[17]+1.623797632095823*fr[13]-1.623797632095823*fl[13]; 
  incr1[18] = (-1.5*fr[18])+1.5*fl[18]+2.662676050517599*fr[10]+2.662676050517599*fl[10]-2.096313728906053*fr[5]+2.096313728906053*fl[5]; 
  incr1[19] = (-2.0625*fr[19])-2.0625*fl[19]+1.623797632095823*fr[15]-1.623797632095823*fl[15]; 

  incr2[2] = 0.4236075534914363*fr[8]+0.4236075534914363*fl[8]-0.609375*fr[2]+0.609375*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[4] = 0.4236075534914363*fr[12]+0.4236075534914363*fl[12]-0.609375*fr[4]+0.609375*fl[4]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = 0.4236075534914363*fr[14]+0.4236075534914363*fl[14]-0.609375*fr[6]+0.609375*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-1.640625*fr[8])-1.640625*fl[8]+2.360099226595144*fr[2]-2.360099226595144*fl[2]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 
  incr2[10] = 0.4236075534914363*fr[18]+0.4236075534914363*fl[18]-0.609375*fr[10]+0.609375*fl[10]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[11] = (-0.609375*fr[11])+0.609375*fl[11]+0.4330127018922194*fr[7]+0.4330127018922194*fl[7]; 
  incr2[12] = (-1.640625*fr[12])-1.640625*fl[12]+2.360099226595145*fr[4]-2.360099226595145*fl[4]-1.677050983124842*fr[1]-1.677050983124842*fl[1]; 
  incr2[14] = (-1.640625*fr[14])-1.640625*fl[14]+2.360099226595145*fr[6]-2.360099226595145*fl[6]-1.677050983124842*fr[3]-1.677050983124842*fl[3]; 
  incr2[16] = (-0.609375*fr[16])+0.609375*fl[16]+0.4330127018922194*fr[9]+0.4330127018922194*fl[9]; 
  incr2[17] = (-0.609375*fr[17])+0.609375*fl[17]+0.4330127018922194*fr[13]+0.4330127018922194*fl[13]; 
  incr2[18] = (-1.640625*fr[18])-1.640625*fl[18]+2.360099226595144*fr[10]-2.360099226595144*fl[10]-1.677050983124842*fr[5]-1.677050983124842*fl[5]; 
  incr2[19] = (-0.609375*fr[19])+0.609375*fl[19]+0.4330127018922194*fr[15]+0.4330127018922194*fl[15]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr1[15]*rdxSq4nur; 
  outr[16] += incr2[16]*rdxSq4nur+incr1[16]*rdxSq4nur; 
  outr[17] += incr2[17]*rdxSq4nur+incr1[17]*rdxSq4nur; 
  outr[18] += incr2[18]*rdxSq4nur+incr1[18]*rdxSq4nur; 
  outr[19] += incr2[19]*rdxSq4nur+incr1[19]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr2[8]*rdxSq4nul-1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr2[12]*rdxSq4nul-1.0*incr1[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr2[14]*rdxSq4nul-1.0*incr1[14]*rdxSq4nul; 
  outl[15] += -1.0*incr1[15]*rdxSq4nul; 
  outl[16] += incr1[16]*rdxSq4nul-1.0*incr2[16]*rdxSq4nul; 
  outl[17] += incr1[17]*rdxSq4nul-1.0*incr2[17]*rdxSq4nul; 
  outl[18] += incr2[18]*rdxSq4nul-1.0*incr1[18]*rdxSq4nul; 
  outl[19] += incr1[19]*rdxSq4nul-1.0*incr2[19]*rdxSq4nul; 

} 
void ConstDiffusionSurf3xSer_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxvr[2]*dxvr[2]); 

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
void ConstDiffusionSurf3xSer_Z_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dxv[3]:    Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxvl[2]*dxvl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxvr[2]*dxvr[2]); 

  double incr1[20]; 
  double incr2[20]; 

  incr1[0] = (-0.6708203932499369*fr[9])+0.6708203932499369*fl[9]+1.190784930203603*fr[3]+1.190784930203603*fl[3]-0.9375*fr[0]+0.9375*fl[0]; 
  incr1[1] = (-0.6708203932499369*fr[15])+0.6708203932499369*fl[15]+1.190784930203603*fr[5]+1.190784930203603*fl[5]-0.9375*fr[1]+0.9375*fl[1]; 
  incr1[2] = (-0.6708203932499369*fr[16])+0.6708203932499369*fl[16]+1.190784930203603*fr[6]+1.190784930203603*fl[6]-0.9375*fr[2]+0.9375*fl[2]; 
  incr1[3] = 1.161895003862225*fr[9]-1.161895003862225*fl[9]-2.0625*fr[3]-2.0625*fl[3]+1.623797632095822*fr[0]-1.623797632095822*fl[0]; 
  incr1[4] = (-0.6708203932499369*fr[19])+0.6708203932499369*fl[19]+1.190784930203603*fr[10]+1.190784930203603*fl[10]-0.9375*fr[4]+0.9375*fl[4]; 
  incr1[5] = 1.161895003862225*fr[15]-1.161895003862225*fl[15]-2.0625*fr[5]-2.0625*fl[5]+1.623797632095822*fr[1]-1.623797632095822*fl[1]; 
  incr1[6] = 1.161895003862225*fr[16]-1.161895003862225*fl[16]-2.0625*fr[6]-2.0625*fl[6]+1.623797632095822*fr[2]-1.623797632095822*fl[2]; 
  incr1[7] = 1.190784930203603*fr[13]+1.190784930203603*fl[13]-0.9375*fr[7]+0.9375*fl[7]; 
  incr1[8] = 1.190784930203603*fr[14]+1.190784930203603*fl[14]-0.9375*fr[8]+0.9375*fl[8]; 
  incr1[9] = (-1.5*fr[9])+1.5*fl[9]+2.662676050517599*fr[3]+2.662676050517599*fl[3]-2.096313728906053*fr[0]+2.096313728906053*fl[0]; 
  incr1[10] = 1.161895003862225*fr[19]-1.161895003862225*fl[19]-2.0625*fr[10]-2.0625*fl[10]+1.623797632095822*fr[4]-1.623797632095822*fl[4]; 
  incr1[11] = 1.190784930203603*fr[17]+1.190784930203603*fl[17]-0.9375*fr[11]+0.9375*fl[11]; 
  incr1[12] = 1.190784930203603*fr[18]+1.190784930203603*fl[18]-0.9375*fr[12]+0.9375*fl[12]; 
  incr1[13] = (-2.0625*fr[13])-2.0625*fl[13]+1.623797632095823*fr[7]-1.623797632095823*fl[7]; 
  incr1[14] = (-2.0625*fr[14])-2.0625*fl[14]+1.623797632095823*fr[8]-1.623797632095823*fl[8]; 
  incr1[15] = (-1.5*fr[15])+1.5*fl[15]+2.662676050517599*fr[5]+2.662676050517599*fl[5]-2.096313728906053*fr[1]+2.096313728906053*fl[1]; 
  incr1[16] = (-1.5*fr[16])+1.5*fl[16]+2.662676050517599*fr[6]+2.662676050517599*fl[6]-2.096313728906053*fr[2]+2.096313728906053*fl[2]; 
  incr1[17] = (-2.0625*fr[17])-2.0625*fl[17]+1.623797632095823*fr[11]-1.623797632095823*fl[11]; 
  incr1[18] = (-2.0625*fr[18])-2.0625*fl[18]+1.623797632095823*fr[12]-1.623797632095823*fl[12]; 
  incr1[19] = (-1.5*fr[19])+1.5*fl[19]+2.662676050517599*fr[10]+2.662676050517599*fl[10]-2.096313728906053*fr[4]+2.096313728906053*fl[4]; 

  incr2[3] = 0.4236075534914363*fr[9]+0.4236075534914363*fl[9]-0.609375*fr[3]+0.609375*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = 0.4236075534914363*fr[15]+0.4236075534914363*fl[15]-0.609375*fr[5]+0.609375*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[6] = 0.4236075534914363*fr[16]+0.4236075534914363*fl[16]-0.609375*fr[6]+0.609375*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[9] = (-1.640625*fr[9])-1.640625*fl[9]+2.360099226595144*fr[3]-2.360099226595144*fl[3]-1.677050983124842*fr[0]-1.677050983124842*fl[0]; 
  incr2[10] = 0.4236075534914363*fr[19]+0.4236075534914363*fl[19]-0.609375*fr[10]+0.609375*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[13] = (-0.609375*fr[13])+0.609375*fl[13]+0.4330127018922194*fr[7]+0.4330127018922194*fl[7]; 
  incr2[14] = (-0.609375*fr[14])+0.609375*fl[14]+0.4330127018922194*fr[8]+0.4330127018922194*fl[8]; 
  incr2[15] = (-1.640625*fr[15])-1.640625*fl[15]+2.360099226595145*fr[5]-2.360099226595145*fl[5]-1.677050983124842*fr[1]-1.677050983124842*fl[1]; 
  incr2[16] = (-1.640625*fr[16])-1.640625*fl[16]+2.360099226595145*fr[6]-2.360099226595145*fl[6]-1.677050983124842*fr[2]-1.677050983124842*fl[2]; 
  incr2[17] = (-0.609375*fr[17])+0.609375*fl[17]+0.4330127018922194*fr[11]+0.4330127018922194*fl[11]; 
  incr2[18] = (-0.609375*fr[18])+0.609375*fl[18]+0.4330127018922194*fr[12]+0.4330127018922194*fl[12]; 
  incr2[19] = (-1.640625*fr[19])-1.640625*fl[19]+2.360099226595144*fr[10]-2.360099226595144*fl[10]-1.677050983124842*fr[4]-1.677050983124842*fl[4]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 
  outr[16] += incr2[16]*rdxSq4nur+incr1[16]*rdxSq4nur; 
  outr[17] += incr2[17]*rdxSq4nur+incr1[17]*rdxSq4nur; 
  outr[18] += incr2[18]*rdxSq4nur+incr1[18]*rdxSq4nur; 
  outr[19] += incr2[19]*rdxSq4nur+incr1[19]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr2[9]*rdxSq4nul-1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr2[15]*rdxSq4nul-1.0*incr1[15]*rdxSq4nul; 
  outl[16] += incr2[16]*rdxSq4nul-1.0*incr1[16]*rdxSq4nul; 
  outl[17] += incr1[17]*rdxSq4nul-1.0*incr2[17]*rdxSq4nul; 
  outl[18] += incr1[18]*rdxSq4nul-1.0*incr2[18]*rdxSq4nul; 
  outl[19] += incr2[19]*rdxSq4nul-1.0*incr1[19]*rdxSq4nul; 

} 
