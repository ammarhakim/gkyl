#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf3xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

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
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 

} 
void ConstDiffusionSurf3xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

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
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += -1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 

} 
void ConstDiffusionSurf3xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

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
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 

  incr1[0] = 6.708203932499369*fr[7]-6.708203932499369*fl[7]-8.11898816047911*fr[1]-8.11898816047911*fl[1]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = (-11.61895003862225*fr[7])+11.61895003862225*fl[7]+14.0625*fr[1]+14.0625*fl[1]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[2] = 6.708203932499369*fr[11]-6.708203932499369*fl[11]-8.11898816047911*fr[4]-8.11898816047911*fl[4]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = 6.708203932499369*fr[13]-6.708203932499369*fl[13]-8.11898816047911*fr[5]-8.11898816047911*fl[5]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = (-11.61895003862225*fr[11])+11.61895003862225*fl[11]+14.0625*fr[4]+14.0625*fl[4]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[5] = (-11.61895003862225*fr[13])+11.61895003862225*fl[13]+14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[6] = 6.708203932499369*fr[17]-6.708203932499369*fl[17]-8.11898816047911*fr[10]-8.11898816047911*fl[10]+4.6875*fr[6]-4.6875*fl[6]; 
  incr1[7] = 15.0*fr[7]-15.0*fl[7]-18.15460943534727*fr[1]-18.15460943534727*fl[1]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[8] = (-8.118988160479114*fr[12])-8.118988160479114*fl[12]+4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = (-8.118988160479114*fr[15])-8.118988160479114*fl[15]+4.6875*fr[9]-4.6875*fl[9]; 
  incr1[10] = (-11.61895003862225*fr[17])+11.61895003862225*fl[17]+14.0625*fr[10]+14.0625*fl[10]-8.11898816047911*fr[6]+8.11898816047911*fl[6]; 
  incr1[11] = 15.0*fr[11]-15.0*fl[11]-18.15460943534727*fr[4]-18.15460943534727*fl[4]+10.48156864453026*fr[2]-10.48156864453026*fl[2]; 
  incr1[12] = 14.0625*fr[12]+14.0625*fl[12]-8.118988160479114*fr[8]+8.118988160479114*fl[8]; 
  incr1[13] = 15.0*fr[13]-15.0*fl[13]-18.15460943534727*fr[5]-18.15460943534727*fl[5]+10.48156864453026*fr[3]-10.48156864453026*fl[3]; 
  incr1[14] = (-8.118988160479114*fr[18])-8.118988160479114*fl[18]+4.6875*fr[14]-4.6875*fl[14]; 
  incr1[15] = 14.0625*fr[15]+14.0625*fl[15]-8.118988160479114*fr[9]+8.118988160479114*fl[9]; 
  incr1[16] = (-8.118988160479114*fr[19])-8.118988160479114*fl[19]+4.6875*fr[16]-4.6875*fl[16]; 
  incr1[17] = 15.0*fr[17]-15.0*fl[17]-18.15460943534727*fr[10]-18.15460943534727*fl[10]+10.48156864453027*fr[6]-10.48156864453027*fl[6]; 
  incr1[18] = 14.0625*fr[18]+14.0625*fl[18]-8.118988160479114*fr[14]+8.118988160479114*fl[14]; 
  incr1[19] = 14.0625*fr[19]+14.0625*fl[19]-8.118988160479114*fr[16]+8.118988160479114*fl[16]; 

  incr2[1] = (-2.541645320948617*fr[7])-2.541645320948617*fl[7]+1.40625*fr[1]-1.40625*fl[1]; 
  incr2[4] = (-2.541645320948617*fr[11])-2.541645320948617*fl[11]+1.40625*fr[4]-1.40625*fl[4]; 
  incr2[5] = (-2.541645320948617*fr[13])-2.541645320948617*fl[13]+1.40625*fr[5]-1.40625*fl[5]; 
  incr2[7] = 9.84375*fr[7]+9.84375*fl[7]-5.446382830604179*fr[1]+5.446382830604179*fl[1]; 
  incr2[10] = (-2.541645320948617*fr[17])-2.541645320948617*fl[17]+1.40625*fr[10]-1.40625*fl[10]; 
  incr2[11] = 9.84375*fr[11]+9.84375*fl[11]-5.446382830604181*fr[4]+5.446382830604181*fl[4]; 
  incr2[12] = 1.40625*fr[12]-1.40625*fl[12]; 
  incr2[13] = 9.84375*fr[13]+9.84375*fl[13]-5.446382830604181*fr[5]+5.446382830604181*fl[5]; 
  incr2[15] = 1.40625*fr[15]-1.40625*fl[15]; 
  incr2[17] = 9.84375*fr[17]+9.84375*fl[17]-5.446382830604179*fr[10]+5.446382830604179*fl[10]; 
  incr2[18] = 1.40625*fr[18]-1.40625*fl[18]; 
  incr2[19] = 1.40625*fr[19]-1.40625*fl[19]; 

  incr3[7] = (-4.5*fr[7])+4.5*fl[7]+7.988028151552797*fr[1]+7.988028151552797*fl[1]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 
  incr3[11] = (-4.5*fr[11])+4.5*fl[11]+7.988028151552798*fr[4]+7.988028151552798*fl[4]-6.288941186718158*fr[2]+6.288941186718158*fl[2]; 
  incr3[13] = (-4.5*fr[13])+4.5*fl[13]+7.988028151552798*fr[5]+7.988028151552798*fl[5]-6.288941186718158*fr[3]+6.288941186718158*fl[3]; 
  incr3[17] = (-4.5*fr[17])+4.5*fl[17]+7.988028151552797*fr[10]+7.988028151552797*fl[10]-6.288941186718159*fr[6]+6.288941186718159*fl[6]; 


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
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += (-1.0*incr3[11]*rdxFnur)-1.0*incr2[11]*rdxFnur-1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur-1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 
  outr[16] += -1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr3[17]*rdxFnur)-1.0*incr2[17]*rdxFnur-1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr2[18]*rdxFnur)-1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr2[19]*rdxFnur)-1.0*incr1[19]*rdxFnur; 

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
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr2[11]*rdxFnul+incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul+incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul; 
  outl[17] += incr3[17]*rdxFnul-1.0*incr2[17]*rdxFnul+incr1[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 

  incr1[0] = 6.708203932499369*fr[8]-6.708203932499369*fl[8]-8.11898816047911*fr[2]-8.11898816047911*fl[2]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = 6.708203932499369*fr[12]-6.708203932499369*fl[12]-8.11898816047911*fr[4]-8.11898816047911*fl[4]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = (-11.61895003862225*fr[8])+11.61895003862225*fl[8]+14.0625*fr[2]+14.0625*fl[2]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[3] = 6.708203932499369*fr[14]-6.708203932499369*fl[14]-8.11898816047911*fr[6]-8.11898816047911*fl[6]+4.6875*fr[3]-4.6875*fl[3]; 
  incr1[4] = (-11.61895003862225*fr[12])+11.61895003862225*fl[12]+14.0625*fr[4]+14.0625*fl[4]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[5] = 6.708203932499369*fr[18]-6.708203932499369*fl[18]-8.11898816047911*fr[10]-8.11898816047911*fl[10]+4.6875*fr[5]-4.6875*fl[5]; 
  incr1[6] = (-11.61895003862225*fr[14])+11.61895003862225*fl[14]+14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[3]+8.11898816047911*fl[3]; 
  incr1[7] = (-8.118988160479114*fr[11])-8.118988160479114*fl[11]+4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = 15.0*fr[8]-15.0*fl[8]-18.15460943534727*fr[2]-18.15460943534727*fl[2]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[9] = (-8.118988160479114*fr[16])-8.118988160479114*fl[16]+4.6875*fr[9]-4.6875*fl[9]; 
  incr1[10] = (-11.61895003862225*fr[18])+11.61895003862225*fl[18]+14.0625*fr[10]+14.0625*fl[10]-8.11898816047911*fr[5]+8.11898816047911*fl[5]; 
  incr1[11] = 14.0625*fr[11]+14.0625*fl[11]-8.118988160479114*fr[7]+8.118988160479114*fl[7]; 
  incr1[12] = 15.0*fr[12]-15.0*fl[12]-18.15460943534727*fr[4]-18.15460943534727*fl[4]+10.48156864453026*fr[1]-10.48156864453026*fl[1]; 
  incr1[13] = (-8.118988160479114*fr[17])-8.118988160479114*fl[17]+4.6875*fr[13]-4.6875*fl[13]; 
  incr1[14] = 15.0*fr[14]-15.0*fl[14]-18.15460943534727*fr[6]-18.15460943534727*fl[6]+10.48156864453026*fr[3]-10.48156864453026*fl[3]; 
  incr1[15] = (-8.118988160479114*fr[19])-8.118988160479114*fl[19]+4.6875*fr[15]-4.6875*fl[15]; 
  incr1[16] = 14.0625*fr[16]+14.0625*fl[16]-8.118988160479114*fr[9]+8.118988160479114*fl[9]; 
  incr1[17] = 14.0625*fr[17]+14.0625*fl[17]-8.118988160479114*fr[13]+8.118988160479114*fl[13]; 
  incr1[18] = 15.0*fr[18]-15.0*fl[18]-18.15460943534727*fr[10]-18.15460943534727*fl[10]+10.48156864453027*fr[5]-10.48156864453027*fl[5]; 
  incr1[19] = 14.0625*fr[19]+14.0625*fl[19]-8.118988160479114*fr[15]+8.118988160479114*fl[15]; 

  incr2[2] = (-2.541645320948617*fr[8])-2.541645320948617*fl[8]+1.40625*fr[2]-1.40625*fl[2]; 
  incr2[4] = (-2.541645320948617*fr[12])-2.541645320948617*fl[12]+1.40625*fr[4]-1.40625*fl[4]; 
  incr2[6] = (-2.541645320948617*fr[14])-2.541645320948617*fl[14]+1.40625*fr[6]-1.40625*fl[6]; 
  incr2[8] = 9.84375*fr[8]+9.84375*fl[8]-5.446382830604179*fr[2]+5.446382830604179*fl[2]; 
  incr2[10] = (-2.541645320948617*fr[18])-2.541645320948617*fl[18]+1.40625*fr[10]-1.40625*fl[10]; 
  incr2[11] = 1.40625*fr[11]-1.40625*fl[11]; 
  incr2[12] = 9.84375*fr[12]+9.84375*fl[12]-5.446382830604181*fr[4]+5.446382830604181*fl[4]; 
  incr2[14] = 9.84375*fr[14]+9.84375*fl[14]-5.446382830604181*fr[6]+5.446382830604181*fl[6]; 
  incr2[16] = 1.40625*fr[16]-1.40625*fl[16]; 
  incr2[17] = 1.40625*fr[17]-1.40625*fl[17]; 
  incr2[18] = 9.84375*fr[18]+9.84375*fl[18]-5.446382830604179*fr[10]+5.446382830604179*fl[10]; 
  incr2[19] = 1.40625*fr[19]-1.40625*fl[19]; 

  incr3[8] = (-4.5*fr[8])+4.5*fl[8]+7.988028151552797*fr[2]+7.988028151552797*fl[2]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 
  incr3[12] = (-4.5*fr[12])+4.5*fl[12]+7.988028151552798*fr[4]+7.988028151552798*fl[4]-6.288941186718158*fr[1]+6.288941186718158*fl[1]; 
  incr3[14] = (-4.5*fr[14])+4.5*fl[14]+7.988028151552798*fr[6]+7.988028151552798*fl[6]-6.288941186718158*fr[3]+6.288941186718158*fl[3]; 
  incr3[18] = (-4.5*fr[18])+4.5*fl[18]+7.988028151552797*fr[10]+7.988028151552797*fl[10]-6.288941186718159*fr[5]+6.288941186718159*fl[5]; 


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
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr3[12]*rdxFnur)-1.0*incr2[12]*rdxFnur-1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur-1.0*incr1[14]*rdxFnur; 
  outr[15] += -1.0*incr1[15]*rdxFnur; 
  outr[16] += (-1.0*incr2[16]*rdxFnur)-1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr2[17]*rdxFnur)-1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr3[18]*rdxFnur)-1.0*incr2[18]*rdxFnur-1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr2[19]*rdxFnur)-1.0*incr1[19]*rdxFnur; 

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
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr2[12]*rdxFnul+incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul+incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr3[18]*rdxFnul-1.0*incr2[18]*rdxFnul+incr1[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf3xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 

  incr1[0] = 6.708203932499369*fr[9]-6.708203932499369*fl[9]-8.11898816047911*fr[3]-8.11898816047911*fl[3]+4.6875*fr[0]-4.6875*fl[0]; 
  incr1[1] = 6.708203932499369*fr[15]-6.708203932499369*fl[15]-8.11898816047911*fr[5]-8.11898816047911*fl[5]+4.6875*fr[1]-4.6875*fl[1]; 
  incr1[2] = 6.708203932499369*fr[16]-6.708203932499369*fl[16]-8.11898816047911*fr[6]-8.11898816047911*fl[6]+4.6875*fr[2]-4.6875*fl[2]; 
  incr1[3] = (-11.61895003862225*fr[9])+11.61895003862225*fl[9]+14.0625*fr[3]+14.0625*fl[3]-8.11898816047911*fr[0]+8.11898816047911*fl[0]; 
  incr1[4] = 6.708203932499369*fr[19]-6.708203932499369*fl[19]-8.11898816047911*fr[10]-8.11898816047911*fl[10]+4.6875*fr[4]-4.6875*fl[4]; 
  incr1[5] = (-11.61895003862225*fr[15])+11.61895003862225*fl[15]+14.0625*fr[5]+14.0625*fl[5]-8.11898816047911*fr[1]+8.11898816047911*fl[1]; 
  incr1[6] = (-11.61895003862225*fr[16])+11.61895003862225*fl[16]+14.0625*fr[6]+14.0625*fl[6]-8.11898816047911*fr[2]+8.11898816047911*fl[2]; 
  incr1[7] = (-8.118988160479114*fr[13])-8.118988160479114*fl[13]+4.6875*fr[7]-4.6875*fl[7]; 
  incr1[8] = (-8.118988160479114*fr[14])-8.118988160479114*fl[14]+4.6875*fr[8]-4.6875*fl[8]; 
  incr1[9] = 15.0*fr[9]-15.0*fl[9]-18.15460943534727*fr[3]-18.15460943534727*fl[3]+10.48156864453027*fr[0]-10.48156864453027*fl[0]; 
  incr1[10] = (-11.61895003862225*fr[19])+11.61895003862225*fl[19]+14.0625*fr[10]+14.0625*fl[10]-8.11898816047911*fr[4]+8.11898816047911*fl[4]; 
  incr1[11] = (-8.118988160479114*fr[17])-8.118988160479114*fl[17]+4.6875*fr[11]-4.6875*fl[11]; 
  incr1[12] = (-8.118988160479114*fr[18])-8.118988160479114*fl[18]+4.6875*fr[12]-4.6875*fl[12]; 
  incr1[13] = 14.0625*fr[13]+14.0625*fl[13]-8.118988160479114*fr[7]+8.118988160479114*fl[7]; 
  incr1[14] = 14.0625*fr[14]+14.0625*fl[14]-8.118988160479114*fr[8]+8.118988160479114*fl[8]; 
  incr1[15] = 15.0*fr[15]-15.0*fl[15]-18.15460943534727*fr[5]-18.15460943534727*fl[5]+10.48156864453026*fr[1]-10.48156864453026*fl[1]; 
  incr1[16] = 15.0*fr[16]-15.0*fl[16]-18.15460943534727*fr[6]-18.15460943534727*fl[6]+10.48156864453026*fr[2]-10.48156864453026*fl[2]; 
  incr1[17] = 14.0625*fr[17]+14.0625*fl[17]-8.118988160479114*fr[11]+8.118988160479114*fl[11]; 
  incr1[18] = 14.0625*fr[18]+14.0625*fl[18]-8.118988160479114*fr[12]+8.118988160479114*fl[12]; 
  incr1[19] = 15.0*fr[19]-15.0*fl[19]-18.15460943534727*fr[10]-18.15460943534727*fl[10]+10.48156864453027*fr[4]-10.48156864453027*fl[4]; 

  incr2[3] = (-2.541645320948617*fr[9])-2.541645320948617*fl[9]+1.40625*fr[3]-1.40625*fl[3]; 
  incr2[5] = (-2.541645320948617*fr[15])-2.541645320948617*fl[15]+1.40625*fr[5]-1.40625*fl[5]; 
  incr2[6] = (-2.541645320948617*fr[16])-2.541645320948617*fl[16]+1.40625*fr[6]-1.40625*fl[6]; 
  incr2[9] = 9.84375*fr[9]+9.84375*fl[9]-5.446382830604179*fr[3]+5.446382830604179*fl[3]; 
  incr2[10] = (-2.541645320948617*fr[19])-2.541645320948617*fl[19]+1.40625*fr[10]-1.40625*fl[10]; 
  incr2[13] = 1.40625*fr[13]-1.40625*fl[13]; 
  incr2[14] = 1.40625*fr[14]-1.40625*fl[14]; 
  incr2[15] = 9.84375*fr[15]+9.84375*fl[15]-5.446382830604181*fr[5]+5.446382830604181*fl[5]; 
  incr2[16] = 9.84375*fr[16]+9.84375*fl[16]-5.446382830604181*fr[6]+5.446382830604181*fl[6]; 
  incr2[17] = 1.40625*fr[17]-1.40625*fl[17]; 
  incr2[18] = 1.40625*fr[18]-1.40625*fl[18]; 
  incr2[19] = 9.84375*fr[19]+9.84375*fl[19]-5.446382830604179*fr[10]+5.446382830604179*fl[10]; 

  incr3[9] = (-4.5*fr[9])+4.5*fl[9]+7.988028151552797*fr[3]+7.988028151552797*fl[3]-6.288941186718159*fr[0]+6.288941186718159*fl[0]; 
  incr3[15] = (-4.5*fr[15])+4.5*fl[15]+7.988028151552798*fr[5]+7.988028151552798*fl[5]-6.288941186718158*fr[1]+6.288941186718158*fl[1]; 
  incr3[16] = (-4.5*fr[16])+4.5*fl[16]+7.988028151552798*fr[6]+7.988028151552798*fl[6]-6.288941186718158*fr[2]+6.288941186718158*fl[2]; 
  incr3[19] = (-4.5*fr[19])+4.5*fl[19]+7.988028151552797*fr[10]+7.988028151552797*fl[10]-6.288941186718159*fr[4]+6.288941186718159*fl[4]; 


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
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr3[15]*rdxFnur)-1.0*incr2[15]*rdxFnur-1.0*incr1[15]*rdxFnur; 
  outr[16] += (-1.0*incr3[16]*rdxFnur)-1.0*incr2[16]*rdxFnur-1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr2[17]*rdxFnur)-1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr2[18]*rdxFnur)-1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr3[19]*rdxFnur)-1.0*incr2[19]*rdxFnur-1.0*incr1[19]*rdxFnur; 

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
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr3[15]*rdxFnul-1.0*incr2[15]*rdxFnul+incr1[15]*rdxFnul; 
  outl[16] += incr3[16]*rdxFnul-1.0*incr2[16]*rdxFnul+incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr3[19]*rdxFnul-1.0*incr2[19]*rdxFnul+incr1[19]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 
  double incr5[20]; 
  double incr6[20]; 

  incr1[0] = (-35.21807064562169*fr[7])+35.21807064562169*fl[7]+34.09975027401226*fr[1]+34.09975027401226*fl[1]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = 60.9994877027668*fr[7]-60.9994877027668*fl[7]-59.0625*fr[1]-59.0625*fl[1]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[2] = (-35.21807064562168*fr[11])+35.21807064562168*fl[11]+34.09975027401226*fr[4]+34.09975027401226*fl[4]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = (-35.21807064562168*fr[13])+35.21807064562168*fl[13]+34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = 60.99948770276682*fr[11]-60.99948770276682*fl[11]-59.0625*fr[4]-59.0625*fl[4]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[5] = 60.99948770276682*fr[13]-60.99948770276682*fl[13]-59.0625*fr[5]-59.0625*fl[5]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[6] = (-35.21807064562169*fr[17])+35.21807064562169*fl[17]+34.09975027401226*fr[10]+34.09975027401226*fl[10]-19.6875*fr[6]+19.6875*fl[6]; 
  incr1[7] = (-78.75*fr[7])+78.75*fl[7]+76.2493596284585*fr[1]+76.2493596284585*fl[1]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[8] = 34.09975027401227*fr[12]+34.09975027401227*fl[12]-19.6875*fr[8]+19.6875*fl[8]; 
  incr1[9] = 34.09975027401227*fr[15]+34.09975027401227*fl[15]-19.6875*fr[9]+19.6875*fl[9]; 
  incr1[10] = 60.9994877027668*fr[17]-60.9994877027668*fl[17]-59.0625*fr[10]-59.0625*fl[10]+34.09975027401226*fr[6]-34.09975027401226*fl[6]; 
  incr1[11] = (-78.75*fr[11])+78.75*fl[11]+76.24935962845854*fr[4]+76.24935962845854*fl[4]-44.02258830702711*fr[2]+44.02258830702711*fl[2]; 
  incr1[12] = (-59.0625*fr[12])-59.0625*fl[12]+34.09975027401227*fr[8]-34.09975027401227*fl[8]; 
  incr1[13] = (-78.75*fr[13])+78.75*fl[13]+76.24935962845854*fr[5]+76.24935962845854*fl[5]-44.02258830702711*fr[3]+44.02258830702711*fl[3]; 
  incr1[14] = 34.09975027401227*fr[18]+34.09975027401227*fl[18]-19.6875*fr[14]+19.6875*fl[14]; 
  incr1[15] = (-59.0625*fr[15])-59.0625*fl[15]+34.09975027401227*fr[9]-34.09975027401227*fl[9]; 
  incr1[16] = 34.09975027401227*fr[19]+34.09975027401227*fl[19]-19.6875*fr[16]+19.6875*fl[16]; 
  incr1[17] = (-78.75*fr[17])+78.75*fl[17]+76.2493596284585*fr[10]+76.2493596284585*fl[10]-44.02258830702712*fr[6]+44.02258830702712*fl[6]; 
  incr1[18] = (-59.0625*fr[18])-59.0625*fl[18]+34.09975027401227*fr[14]-34.09975027401227*fl[14]; 
  incr1[19] = (-59.0625*fr[19])-59.0625*fl[19]+34.09975027401227*fr[16]-34.09975027401227*fl[16]; 

  incr2[1] = 9.531169953557313*fr[7]+9.531169953557313*fl[7]-2.4609375*fr[1]+2.4609375*fl[1]; 
  incr2[4] = 9.531169953557317*fr[11]+9.531169953557317*fl[11]-2.4609375*fr[4]+2.4609375*fl[4]; 
  incr2[5] = 9.531169953557317*fr[13]+9.531169953557317*fl[13]-2.4609375*fr[5]+2.4609375*fl[5]; 
  incr2[7] = (-36.9140625*fr[7])-36.9140625*fl[7]+9.531169953557313*fr[1]-9.531169953557313*fl[1]; 
  incr2[10] = 9.531169953557313*fr[17]+9.531169953557313*fl[17]-2.4609375*fr[10]+2.4609375*fl[10]; 
  incr2[11] = (-36.9140625*fr[11])-36.9140625*fl[11]+9.531169953557317*fr[4]-9.531169953557317*fl[4]; 
  incr2[12] = 2.4609375*fl[12]-2.4609375*fr[12]; 
  incr2[13] = (-36.9140625*fr[13])-36.9140625*fl[13]+9.531169953557317*fr[5]-9.531169953557317*fl[5]; 
  incr2[15] = 2.4609375*fl[15]-2.4609375*fr[15]; 
  incr2[17] = (-36.9140625*fr[17])-36.9140625*fl[17]+9.531169953557313*fr[10]-9.531169953557313*fl[10]; 
  incr2[18] = 2.4609375*fl[18]-2.4609375*fr[18]; 
  incr2[19] = 2.4609375*fl[19]-2.4609375*fr[19]; 

  incr3[7] = 45.0*fr[7]-45.0*fl[7]-54.4638283060418*fr[1]-54.4638283060418*fl[1]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 
  incr3[11] = 45.0*fr[11]-45.0*fl[11]-54.46382830604181*fr[4]-54.46382830604181*fl[4]+31.44470593359079*fr[2]-31.44470593359079*fl[2]; 
  incr3[13] = 45.0*fr[13]-45.0*fl[13]-54.46382830604181*fr[5]-54.46382830604181*fl[5]+31.44470593359079*fr[3]-31.44470593359079*fl[3]; 
  incr3[17] = 45.0*fr[17]-45.0*fl[17]-54.4638283060418*fr[10]-54.4638283060418*fl[10]+31.4447059335908*fr[6]-31.4447059335908*fl[6]; 




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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr3[11]*rdxFnur+incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr3[13]*rdxFnur+incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr1[16]*rdxFnur; 
  outr[17] += incr3[17]*rdxFnur+incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += (-1.0*incr3[7]*rdxFnul)+incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += (-1.0*incr3[11]*rdxFnul)+incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += (-1.0*incr3[13]*rdxFnul)+incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr1[16]*rdxFnul; 
  outl[17] += (-1.0*incr3[17]*rdxFnul)+incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 
  double incr5[20]; 
  double incr6[20]; 

  incr1[0] = (-35.21807064562169*fr[8])+35.21807064562169*fl[8]+34.09975027401226*fr[2]+34.09975027401226*fl[2]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = (-35.21807064562168*fr[12])+35.21807064562168*fl[12]+34.09975027401226*fr[4]+34.09975027401226*fl[4]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = 60.9994877027668*fr[8]-60.9994877027668*fl[8]-59.0625*fr[2]-59.0625*fl[2]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[3] = (-35.21807064562168*fr[14])+35.21807064562168*fl[14]+34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[3]+19.6875*fl[3]; 
  incr1[4] = 60.99948770276682*fr[12]-60.99948770276682*fl[12]-59.0625*fr[4]-59.0625*fl[4]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[5] = (-35.21807064562169*fr[18])+35.21807064562169*fl[18]+34.09975027401226*fr[10]+34.09975027401226*fl[10]-19.6875*fr[5]+19.6875*fl[5]; 
  incr1[6] = 60.99948770276682*fr[14]-60.99948770276682*fl[14]-59.0625*fr[6]-59.0625*fl[6]+34.09975027401226*fr[3]-34.09975027401226*fl[3]; 
  incr1[7] = 34.09975027401227*fr[11]+34.09975027401227*fl[11]-19.6875*fr[7]+19.6875*fl[7]; 
  incr1[8] = (-78.75*fr[8])+78.75*fl[8]+76.2493596284585*fr[2]+76.2493596284585*fl[2]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[9] = 34.09975027401227*fr[16]+34.09975027401227*fl[16]-19.6875*fr[9]+19.6875*fl[9]; 
  incr1[10] = 60.9994877027668*fr[18]-60.9994877027668*fl[18]-59.0625*fr[10]-59.0625*fl[10]+34.09975027401226*fr[5]-34.09975027401226*fl[5]; 
  incr1[11] = (-59.0625*fr[11])-59.0625*fl[11]+34.09975027401227*fr[7]-34.09975027401227*fl[7]; 
  incr1[12] = (-78.75*fr[12])+78.75*fl[12]+76.24935962845854*fr[4]+76.24935962845854*fl[4]-44.02258830702711*fr[1]+44.02258830702711*fl[1]; 
  incr1[13] = 34.09975027401227*fr[17]+34.09975027401227*fl[17]-19.6875*fr[13]+19.6875*fl[13]; 
  incr1[14] = (-78.75*fr[14])+78.75*fl[14]+76.24935962845854*fr[6]+76.24935962845854*fl[6]-44.02258830702711*fr[3]+44.02258830702711*fl[3]; 
  incr1[15] = 34.09975027401227*fr[19]+34.09975027401227*fl[19]-19.6875*fr[15]+19.6875*fl[15]; 
  incr1[16] = (-59.0625*fr[16])-59.0625*fl[16]+34.09975027401227*fr[9]-34.09975027401227*fl[9]; 
  incr1[17] = (-59.0625*fr[17])-59.0625*fl[17]+34.09975027401227*fr[13]-34.09975027401227*fl[13]; 
  incr1[18] = (-78.75*fr[18])+78.75*fl[18]+76.2493596284585*fr[10]+76.2493596284585*fl[10]-44.02258830702712*fr[5]+44.02258830702712*fl[5]; 
  incr1[19] = (-59.0625*fr[19])-59.0625*fl[19]+34.09975027401227*fr[15]-34.09975027401227*fl[15]; 

  incr2[2] = 9.531169953557313*fr[8]+9.531169953557313*fl[8]-2.4609375*fr[2]+2.4609375*fl[2]; 
  incr2[4] = 9.531169953557317*fr[12]+9.531169953557317*fl[12]-2.4609375*fr[4]+2.4609375*fl[4]; 
  incr2[6] = 9.531169953557317*fr[14]+9.531169953557317*fl[14]-2.4609375*fr[6]+2.4609375*fl[6]; 
  incr2[8] = (-36.9140625*fr[8])-36.9140625*fl[8]+9.531169953557313*fr[2]-9.531169953557313*fl[2]; 
  incr2[10] = 9.531169953557313*fr[18]+9.531169953557313*fl[18]-2.4609375*fr[10]+2.4609375*fl[10]; 
  incr2[11] = 2.4609375*fl[11]-2.4609375*fr[11]; 
  incr2[12] = (-36.9140625*fr[12])-36.9140625*fl[12]+9.531169953557317*fr[4]-9.531169953557317*fl[4]; 
  incr2[14] = (-36.9140625*fr[14])-36.9140625*fl[14]+9.531169953557317*fr[6]-9.531169953557317*fl[6]; 
  incr2[16] = 2.4609375*fl[16]-2.4609375*fr[16]; 
  incr2[17] = 2.4609375*fl[17]-2.4609375*fr[17]; 
  incr2[18] = (-36.9140625*fr[18])-36.9140625*fl[18]+9.531169953557313*fr[10]-9.531169953557313*fl[10]; 
  incr2[19] = 2.4609375*fl[19]-2.4609375*fr[19]; 

  incr3[8] = 45.0*fr[8]-45.0*fl[8]-54.4638283060418*fr[2]-54.4638283060418*fl[2]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 
  incr3[12] = 45.0*fr[12]-45.0*fl[12]-54.46382830604181*fr[4]-54.46382830604181*fl[4]+31.44470593359079*fr[1]-31.44470593359079*fl[1]; 
  incr3[14] = 45.0*fr[14]-45.0*fl[14]-54.46382830604181*fr[6]-54.46382830604181*fl[6]+31.44470593359079*fr[3]-31.44470593359079*fl[3]; 
  incr3[18] = 45.0*fr[18]-45.0*fl[18]-54.4638283060418*fr[10]-54.4638283060418*fl[10]+31.4447059335908*fr[5]-31.4447059335908*fl[5]; 




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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr3[12]*rdxFnur+incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr3[14]*rdxFnur+incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr3[18]*rdxFnur+incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += (-1.0*incr3[8]*rdxFnul)+incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += (-1.0*incr3[12]*rdxFnul)+incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += (-1.0*incr3[14]*rdxFnul)+incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += -1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += (-1.0*incr3[18]*rdxFnul)+incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 

} 
void ConstHyperDiffusion6Surf3xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[20]; 
  double incr2[20]; 
  double incr3[20]; 
  double incr4[20]; 
  double incr5[20]; 
  double incr6[20]; 

  incr1[0] = (-35.21807064562169*fr[9])+35.21807064562169*fl[9]+34.09975027401226*fr[3]+34.09975027401226*fl[3]-19.6875*fr[0]+19.6875*fl[0]; 
  incr1[1] = (-35.21807064562168*fr[15])+35.21807064562168*fl[15]+34.09975027401226*fr[5]+34.09975027401226*fl[5]-19.6875*fr[1]+19.6875*fl[1]; 
  incr1[2] = (-35.21807064562168*fr[16])+35.21807064562168*fl[16]+34.09975027401226*fr[6]+34.09975027401226*fl[6]-19.6875*fr[2]+19.6875*fl[2]; 
  incr1[3] = 60.9994877027668*fr[9]-60.9994877027668*fl[9]-59.0625*fr[3]-59.0625*fl[3]+34.09975027401226*fr[0]-34.09975027401226*fl[0]; 
  incr1[4] = (-35.21807064562169*fr[19])+35.21807064562169*fl[19]+34.09975027401226*fr[10]+34.09975027401226*fl[10]-19.6875*fr[4]+19.6875*fl[4]; 
  incr1[5] = 60.99948770276682*fr[15]-60.99948770276682*fl[15]-59.0625*fr[5]-59.0625*fl[5]+34.09975027401226*fr[1]-34.09975027401226*fl[1]; 
  incr1[6] = 60.99948770276682*fr[16]-60.99948770276682*fl[16]-59.0625*fr[6]-59.0625*fl[6]+34.09975027401226*fr[2]-34.09975027401226*fl[2]; 
  incr1[7] = 34.09975027401227*fr[13]+34.09975027401227*fl[13]-19.6875*fr[7]+19.6875*fl[7]; 
  incr1[8] = 34.09975027401227*fr[14]+34.09975027401227*fl[14]-19.6875*fr[8]+19.6875*fl[8]; 
  incr1[9] = (-78.75*fr[9])+78.75*fl[9]+76.2493596284585*fr[3]+76.2493596284585*fl[3]-44.02258830702712*fr[0]+44.02258830702712*fl[0]; 
  incr1[10] = 60.9994877027668*fr[19]-60.9994877027668*fl[19]-59.0625*fr[10]-59.0625*fl[10]+34.09975027401226*fr[4]-34.09975027401226*fl[4]; 
  incr1[11] = 34.09975027401227*fr[17]+34.09975027401227*fl[17]-19.6875*fr[11]+19.6875*fl[11]; 
  incr1[12] = 34.09975027401227*fr[18]+34.09975027401227*fl[18]-19.6875*fr[12]+19.6875*fl[12]; 
  incr1[13] = (-59.0625*fr[13])-59.0625*fl[13]+34.09975027401227*fr[7]-34.09975027401227*fl[7]; 
  incr1[14] = (-59.0625*fr[14])-59.0625*fl[14]+34.09975027401227*fr[8]-34.09975027401227*fl[8]; 
  incr1[15] = (-78.75*fr[15])+78.75*fl[15]+76.24935962845854*fr[5]+76.24935962845854*fl[5]-44.02258830702711*fr[1]+44.02258830702711*fl[1]; 
  incr1[16] = (-78.75*fr[16])+78.75*fl[16]+76.24935962845854*fr[6]+76.24935962845854*fl[6]-44.02258830702711*fr[2]+44.02258830702711*fl[2]; 
  incr1[17] = (-59.0625*fr[17])-59.0625*fl[17]+34.09975027401227*fr[11]-34.09975027401227*fl[11]; 
  incr1[18] = (-59.0625*fr[18])-59.0625*fl[18]+34.09975027401227*fr[12]-34.09975027401227*fl[12]; 
  incr1[19] = (-78.75*fr[19])+78.75*fl[19]+76.2493596284585*fr[10]+76.2493596284585*fl[10]-44.02258830702712*fr[4]+44.02258830702712*fl[4]; 

  incr2[3] = 9.531169953557313*fr[9]+9.531169953557313*fl[9]-2.4609375*fr[3]+2.4609375*fl[3]; 
  incr2[5] = 9.531169953557317*fr[15]+9.531169953557317*fl[15]-2.4609375*fr[5]+2.4609375*fl[5]; 
  incr2[6] = 9.531169953557317*fr[16]+9.531169953557317*fl[16]-2.4609375*fr[6]+2.4609375*fl[6]; 
  incr2[9] = (-36.9140625*fr[9])-36.9140625*fl[9]+9.531169953557313*fr[3]-9.531169953557313*fl[3]; 
  incr2[10] = 9.531169953557313*fr[19]+9.531169953557313*fl[19]-2.4609375*fr[10]+2.4609375*fl[10]; 
  incr2[13] = 2.4609375*fl[13]-2.4609375*fr[13]; 
  incr2[14] = 2.4609375*fl[14]-2.4609375*fr[14]; 
  incr2[15] = (-36.9140625*fr[15])-36.9140625*fl[15]+9.531169953557317*fr[5]-9.531169953557317*fl[5]; 
  incr2[16] = (-36.9140625*fr[16])-36.9140625*fl[16]+9.531169953557317*fr[6]-9.531169953557317*fl[6]; 
  incr2[17] = 2.4609375*fl[17]-2.4609375*fr[17]; 
  incr2[18] = 2.4609375*fl[18]-2.4609375*fr[18]; 
  incr2[19] = (-36.9140625*fr[19])-36.9140625*fl[19]+9.531169953557313*fr[10]-9.531169953557313*fl[10]; 

  incr3[9] = 45.0*fr[9]-45.0*fl[9]-54.4638283060418*fr[3]-54.4638283060418*fl[3]+31.4447059335908*fr[0]-31.4447059335908*fl[0]; 
  incr3[15] = 45.0*fr[15]-45.0*fl[15]-54.46382830604181*fr[5]-54.46382830604181*fl[5]+31.44470593359079*fr[1]-31.44470593359079*fl[1]; 
  incr3[16] = 45.0*fr[16]-45.0*fl[16]-54.46382830604181*fr[6]-54.46382830604181*fl[6]+31.44470593359079*fr[2]-31.44470593359079*fl[2]; 
  incr3[19] = 45.0*fr[19]-45.0*fl[19]-54.4638283060418*fr[10]-54.4638283060418*fl[10]+31.4447059335908*fr[4]-31.4447059335908*fl[4]; 




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
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr3[15]*rdxFnur+incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr3[16]*rdxFnur+incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr3[19]*rdxFnur+incr2[19]*rdxFnur+incr1[19]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += (-1.0*incr3[9]*rdxFnul)+incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += (-1.0*incr3[15]*rdxFnul)+incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 
  outl[16] += (-1.0*incr3[16]*rdxFnul)+incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += (-1.0*incr3[19]*rdxFnul)+incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 

} 
