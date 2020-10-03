#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf4xSerP1_diffDirs2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs23_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs23_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs234_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs234_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs234_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs1234_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs1234_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs1234_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs1234_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs123_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs123_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs123_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs24_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs24_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs124_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs124_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs124_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs12_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs12_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[6] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[13] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr2[2]*rdxSq4nur+incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += incr1[2]*rdxSq4nul-1.0*incr2[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += -1.0*incr1[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs34_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs34_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs134_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs134_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs134_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[2]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[2]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs13_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs13_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[2]*dxl[2]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[8] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[12] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr2[3]*rdxSq4nur+incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr2[7]*rdxSq4nur+incr1[7]*rdxSq4nur; 
  outr[8] += incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += incr1[3]*rdxSq4nul-1.0*incr2[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += incr1[7]*rdxSq4nul-1.0*incr2[7]*rdxSq4nul; 
  outl[8] += -1.0*incr1[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += -1.0*incr1[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs4_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs14_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs14_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[1]/(dxl[3]*dxl[3]); 
  double rdxSq4nur = 4.0*nu[1]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[11] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr2[4]*rdxSq4nur+incr1[4]*rdxSq4nur; 
  outr[5] += incr1[5]*rdxSq4nur; 
  outr[6] += incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr2[9]*rdxSq4nur+incr1[9]*rdxSq4nur; 
  outr[10] += incr2[10]*rdxSq4nur+incr1[10]*rdxSq4nur; 
  outr[11] += incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr2[14]*rdxSq4nur+incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += -1.0*incr1[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += incr1[4]*rdxSq4nul-1.0*incr2[4]*rdxSq4nul; 
  outl[5] += -1.0*incr1[5]*rdxSq4nul; 
  outl[6] += -1.0*incr1[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += incr1[9]*rdxSq4nul-1.0*incr2[9]*rdxSq4nul; 
  outl[10] += incr1[10]*rdxSq4nul-1.0*incr2[10]*rdxSq4nul; 
  outl[11] += -1.0*incr1[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += incr1[14]*rdxSq4nul-1.0*incr2[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
void ConstDiffusionSurf4xSerP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq4nul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxSq4nur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[7] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[9] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[14] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 

  outr[0] += incr1[0]*rdxSq4nur; 
  outr[1] += incr2[1]*rdxSq4nur+incr1[1]*rdxSq4nur; 
  outr[2] += incr1[2]*rdxSq4nur; 
  outr[3] += incr1[3]*rdxSq4nur; 
  outr[4] += incr1[4]*rdxSq4nur; 
  outr[5] += incr2[5]*rdxSq4nur+incr1[5]*rdxSq4nur; 
  outr[6] += incr2[6]*rdxSq4nur+incr1[6]*rdxSq4nur; 
  outr[7] += incr1[7]*rdxSq4nur; 
  outr[8] += incr2[8]*rdxSq4nur+incr1[8]*rdxSq4nur; 
  outr[9] += incr1[9]*rdxSq4nur; 
  outr[10] += incr1[10]*rdxSq4nur; 
  outr[11] += incr2[11]*rdxSq4nur+incr1[11]*rdxSq4nur; 
  outr[12] += incr2[12]*rdxSq4nur+incr1[12]*rdxSq4nur; 
  outr[13] += incr2[13]*rdxSq4nur+incr1[13]*rdxSq4nur; 
  outr[14] += incr1[14]*rdxSq4nur; 
  outr[15] += incr2[15]*rdxSq4nur+incr1[15]*rdxSq4nur; 

  outl[0] += -1.0*incr1[0]*rdxSq4nul; 
  outl[1] += incr1[1]*rdxSq4nul-1.0*incr2[1]*rdxSq4nul; 
  outl[2] += -1.0*incr1[2]*rdxSq4nul; 
  outl[3] += -1.0*incr1[3]*rdxSq4nul; 
  outl[4] += -1.0*incr1[4]*rdxSq4nul; 
  outl[5] += incr1[5]*rdxSq4nul-1.0*incr2[5]*rdxSq4nul; 
  outl[6] += incr1[6]*rdxSq4nul-1.0*incr2[6]*rdxSq4nul; 
  outl[7] += -1.0*incr1[7]*rdxSq4nul; 
  outl[8] += incr1[8]*rdxSq4nul-1.0*incr2[8]*rdxSq4nul; 
  outl[9] += -1.0*incr1[9]*rdxSq4nul; 
  outl[10] += -1.0*incr1[10]*rdxSq4nul; 
  outl[11] += incr1[11]*rdxSq4nul-1.0*incr2[11]*rdxSq4nul; 
  outl[12] += incr1[12]*rdxSq4nul-1.0*incr2[12]*rdxSq4nul; 
  outl[13] += incr1[13]*rdxSq4nul-1.0*incr2[13]*rdxSq4nul; 
  outl[14] += -1.0*incr1[14]*rdxSq4nul; 
  outl[15] += incr1[15]*rdxSq4nul-1.0*incr2[15]*rdxSq4nul; 

} 
