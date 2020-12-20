#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

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
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 

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
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 

} 
void ConstDiffusionSurf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

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
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 

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
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 

} 
void ConstDiffusionSurf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

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
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 

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
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 

} 
void ConstDiffusionSurf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

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
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 

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
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  incr1[0] = (-1.623797632095822*fr[1])-1.623797632095822*fl[1]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = 2.8125*fr[1]+2.8125*fl[1]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[2] = (-1.623797632095822*fr[5])-1.623797632095822*fl[5]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = (-1.623797632095822*fr[6])-1.623797632095822*fl[6]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = (-1.623797632095822*fr[8])-1.623797632095822*fl[8]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = 2.8125*fr[5]+2.8125*fl[5]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[6] = 2.8125*fr[6]+2.8125*fl[6]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[7] = (-1.623797632095822*fr[11])-1.623797632095822*fl[11]+0.9375*fr[7]-0.9375*fl[7]; 
  incr1[8] = 2.8125*fr[8]+2.8125*fl[8]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[9] = (-1.623797632095822*fr[12])-1.623797632095822*fl[12]+0.9375*fr[9]-0.9375*fl[9]; 
  incr1[10] = (-1.623797632095822*fr[13])-1.623797632095822*fl[13]+0.9375*fr[10]-0.9375*fl[10]; 
  incr1[11] = 2.8125*fr[11]+2.8125*fl[11]-1.623797632095822*fr[7]+1.623797632095822*fl[7]; 
  incr1[12] = 2.8125*fr[12]+2.8125*fl[12]-1.623797632095822*fr[9]+1.623797632095822*fl[9]; 
  incr1[13] = 2.8125*fr[13]+2.8125*fl[13]-1.623797632095822*fr[10]+1.623797632095822*fl[10]; 
  incr1[14] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[14]-0.9375*fl[14]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[14]+1.623797632095822*fl[14]; 

  incr2[1] = 0.75*fr[1]-0.75*fl[1]; 
  incr2[5] = 0.75*fr[5]-0.75*fl[5]; 
  incr2[6] = 0.75*fr[6]-0.75*fl[6]; 
  incr2[8] = 0.75*fr[8]-0.75*fl[8]; 
  incr2[11] = 0.75*fr[11]-0.75*fl[11]; 
  incr2[12] = 0.75*fr[12]-0.75*fl[12]; 
  incr2[13] = 0.75*fr[13]-0.75*fl[13]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 



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
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 

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
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  incr1[0] = (-1.623797632095822*fr[2])-1.623797632095822*fl[2]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[5])-1.623797632095822*fl[5]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = 2.8125*fr[2]+2.8125*fl[2]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[3] = (-1.623797632095822*fr[7])-1.623797632095822*fl[7]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = (-1.623797632095822*fr[9])-1.623797632095822*fl[9]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = 2.8125*fr[5]+2.8125*fl[5]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[6] = (-1.623797632095822*fr[11])-1.623797632095822*fl[11]+0.9375*fr[6]-0.9375*fl[6]; 
  incr1[7] = 2.8125*fr[7]+2.8125*fl[7]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[8] = (-1.623797632095822*fr[12])-1.623797632095822*fl[12]+0.9375*fr[8]-0.9375*fl[8]; 
  incr1[9] = 2.8125*fr[9]+2.8125*fl[9]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[10] = (-1.623797632095822*fr[14])-1.623797632095822*fl[14]+0.9375*fr[10]-0.9375*fl[10]; 
  incr1[11] = 2.8125*fr[11]+2.8125*fl[11]-1.623797632095822*fr[6]+1.623797632095822*fl[6]; 
  incr1[12] = 2.8125*fr[12]+2.8125*fl[12]-1.623797632095822*fr[8]+1.623797632095822*fl[8]; 
  incr1[13] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[13]-0.9375*fl[13]; 
  incr1[14] = 2.8125*fr[14]+2.8125*fl[14]-1.623797632095822*fr[10]+1.623797632095822*fl[10]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[13]+1.623797632095822*fl[13]; 

  incr2[2] = 0.75*fr[2]-0.75*fl[2]; 
  incr2[5] = 0.75*fr[5]-0.75*fl[5]; 
  incr2[7] = 0.75*fr[7]-0.75*fl[7]; 
  incr2[9] = 0.75*fr[9]-0.75*fl[9]; 
  incr2[11] = 0.75*fr[11]-0.75*fl[11]; 
  incr2[12] = 0.75*fr[12]-0.75*fl[12]; 
  incr2[14] = 0.75*fr[14]-0.75*fl[14]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 



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
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 

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
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  incr1[0] = (-1.623797632095822*fr[3])-1.623797632095822*fl[3]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[6])-1.623797632095822*fl[6]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = (-1.623797632095822*fr[7])-1.623797632095822*fl[7]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = 2.8125*fr[3]+2.8125*fl[3]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[4] = (-1.623797632095822*fr[10])-1.623797632095822*fl[10]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = (-1.623797632095822*fr[11])-1.623797632095822*fl[11]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = 2.8125*fr[6]+2.8125*fl[6]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[7] = 2.8125*fr[7]+2.8125*fl[7]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[8] = (-1.623797632095822*fr[13])-1.623797632095822*fl[13]+0.9375*fr[8]-0.9375*fl[8]; 
  incr1[9] = (-1.623797632095822*fr[14])-1.623797632095822*fl[14]+0.9375*fr[9]-0.9375*fl[9]; 
  incr1[10] = 2.8125*fr[10]+2.8125*fl[10]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[11] = 2.8125*fr[11]+2.8125*fl[11]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[12] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[12]-0.9375*fl[12]; 
  incr1[13] = 2.8125*fr[13]+2.8125*fl[13]-1.623797632095822*fr[8]+1.623797632095822*fl[8]; 
  incr1[14] = 2.8125*fr[14]+2.8125*fl[14]-1.623797632095822*fr[9]+1.623797632095822*fl[9]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[12]+1.623797632095822*fl[12]; 

  incr2[3] = 0.75*fr[3]-0.75*fl[3]; 
  incr2[6] = 0.75*fr[6]-0.75*fl[6]; 
  incr2[7] = 0.75*fr[7]-0.75*fl[7]; 
  incr2[10] = 0.75*fr[10]-0.75*fl[10]; 
  incr2[11] = 0.75*fr[11]-0.75*fl[11]; 
  incr2[13] = 0.75*fr[13]-0.75*fl[13]; 
  incr2[14] = 0.75*fr[14]-0.75*fl[14]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 



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
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 

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
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  incr1[0] = (-1.623797632095822*fr[4])-1.623797632095822*fl[4]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[8])-1.623797632095822*fl[8]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = (-1.623797632095822*fr[9])-1.623797632095822*fl[9]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = (-1.623797632095822*fr[10])-1.623797632095822*fl[10]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = 2.8125*fr[4]+2.8125*fl[4]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[5] = (-1.623797632095822*fr[12])-1.623797632095822*fl[12]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = (-1.623797632095822*fr[13])-1.623797632095822*fl[13]+0.9375*fr[6]-0.9375*fl[6]; 
  incr1[7] = (-1.623797632095822*fr[14])-1.623797632095822*fl[14]+0.9375*fr[7]-0.9375*fl[7]; 
  incr1[8] = 2.8125*fr[8]+2.8125*fl[8]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[9] = 2.8125*fr[9]+2.8125*fl[9]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[10] = 2.8125*fr[10]+2.8125*fl[10]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[11] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[11]-0.9375*fl[11]; 
  incr1[12] = 2.8125*fr[12]+2.8125*fl[12]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[13] = 2.8125*fr[13]+2.8125*fl[13]-1.623797632095822*fr[6]+1.623797632095822*fl[6]; 
  incr1[14] = 2.8125*fr[14]+2.8125*fl[14]-1.623797632095822*fr[7]+1.623797632095822*fl[7]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[11]+1.623797632095822*fl[11]; 

  incr2[4] = 0.75*fr[4]-0.75*fl[4]; 
  incr2[8] = 0.75*fr[8]-0.75*fl[8]; 
  incr2[9] = 0.75*fr[9]-0.75*fl[9]; 
  incr2[10] = 0.75*fr[10]-0.75*fl[10]; 
  incr2[12] = 0.75*fr[12]-0.75*fl[12]; 
  incr2[13] = 0.75*fr[13]-0.75*fl[13]; 
  incr2[14] = 0.75*fr[14]-0.75*fl[14]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 



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
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 

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
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 

} 
