#include <VlasovModDecl.h> 
double VlasovSurfElcMag2x2vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 

  double incr[5]; 

  double favg[5]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  double fjump[5]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  const double amid = 0.5*abar0[0]; 
  incr[0] = 0.03608439182435161*B2[0]*favg[4]*dv2-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.2165063509461096*abar0[0]*favg[3]+0.125*abar0[2]*favg[2]+0.125*abar0[1]*favg[1]+0.125*abar0[0]*favg[0]; 
  incr[1] = 0.03608439182435161*B2[1]*favg[4]*dv2-0.25*fjump[1]*amax+0.2165063509461096*abar0[1]*favg[3]+0.125*abar0[0]*favg[1]+0.125*favg[0]*abar0[1]; 
  incr[2] = 0.03608439182435161*B2[2]*favg[4]*dv2-0.25*fjump[2]*amax+0.2165063509461096*abar0[2]*favg[3]+0.125*abar0[0]*favg[2]+0.125*favg[0]*abar0[2]; 
  incr[3] = (-0.0625*B2[0]*favg[4]*dv2)+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.375*abar0[0]*favg[3]-0.2165063509461096*abar0[2]*favg[2]-0.2165063509461096*abar0[1]*favg[1]-0.2165063509461096*abar0[0]*favg[0]; 
  incr[4] = 0.0625*B2[0]*favg[3]*dv2+0.03608439182435161*favg[2]*B2[2]*dv2+0.03608439182435161*favg[1]*B2[1]*dv2+0.03608439182435161*favg[0]*B2[0]*dv2-0.25*fjump[4]*amax+0.125*abar0[0]*favg[4]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += -1.0*incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x2vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 
  abar0[3] = E0[3]+wv2*B2[3]; 
  abar0[4] = E0[4]+wv2*B2[4]; 
  abar0[5] = E0[5]+wv2*B2[5]; 

  double incr[15]; 

  double favg[15]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = -1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  double fjump[15]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = -1*fr[6]-fl[6]; 
  fjump[7] = -1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = 1*fr[9]-fl[9]; 
  fjump[10] = -1*fr[10]-fl[10]; 
  fjump[11] = 1*fr[11]-fl[11]; 
  fjump[12] = 1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = 1*fr[14]-fl[14]; 
  const double amid = (-0.5590169943749475*abar0[5])-0.5590169943749475*abar0[4]+0.5*abar0[0]; 
  incr[0] = 0.0625*B2[0]*favg[10]*dv2+0.03608439182435161*B2[2]*favg[9]*dv2+0.03608439182435161*B2[1]*favg[8]*dv2+0.03608439182435161*B2[0]*favg[4]*dv2-0.5590169943749475*fjump[13]*amax-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.2795084971874737*abar0[0]*favg[13]+0.125*abar0[5]*favg[12]+0.125*abar0[4]*favg[11]+0.2165063509461096*abar0[2]*favg[7]+0.2165063509461096*abar0[1]*favg[6]+0.125*abar0[3]*favg[5]+0.2165063509461096*abar0[0]*favg[3]+0.125*abar0[2]*favg[2]+0.125*abar0[1]*favg[1]+0.125*abar0[0]*favg[0]; 
  incr[1] = 0.0625*B2[1]*favg[10]*dv2+0.03608439182435161*B2[3]*favg[9]*dv2+0.03227486121839514*B2[4]*favg[8]*dv2+0.03608439182435161*B2[0]*favg[8]*dv2+0.03608439182435161*B2[1]*favg[4]*dv2-0.4330127018922193*fjump[6]*amax-0.25*fjump[1]*amax+0.2795084971874738*abar0[1]*favg[13]+0.1118033988749895*abar0[1]*favg[11]+0.2165063509461096*abar0[3]*favg[7]+0.1936491673103708*abar0[4]*favg[6]+0.2165063509461096*abar0[0]*favg[6]+0.125*abar0[2]*favg[5]+0.1118033988749895*favg[1]*abar0[4]+0.2165063509461096*abar0[1]*favg[3]+0.125*favg[2]*abar0[3]+0.125*abar0[0]*favg[1]+0.125*favg[0]*abar0[1]; 
  incr[2] = 0.0625*B2[2]*favg[10]*dv2+0.03227486121839514*B2[5]*favg[9]*dv2+0.03608439182435161*B2[0]*favg[9]*dv2+0.03608439182435161*B2[3]*favg[8]*dv2+0.03608439182435161*B2[2]*favg[4]*dv2-0.4330127018922193*fjump[7]*amax-0.25*fjump[2]*amax+0.2795084971874738*abar0[2]*favg[13]+0.1118033988749895*abar0[2]*favg[12]+0.1936491673103708*abar0[5]*favg[7]+0.2165063509461096*abar0[0]*favg[7]+0.2165063509461096*abar0[3]*favg[6]+0.125*abar0[1]*favg[5]+0.1118033988749895*favg[2]*abar0[5]+0.2165063509461096*abar0[2]*favg[3]+0.125*favg[1]*abar0[3]+0.125*abar0[0]*favg[2]+0.125*favg[0]*abar0[2]; 
  incr[3] = (-0.1082531754730548*B2[0]*favg[10]*dv2)-0.0625*B2[2]*favg[9]*dv2-0.0625*B2[1]*favg[8]*dv2-0.0625*B2[0]*favg[4]*dv2+0.9682458365518543*fjump[13]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.4841229182759271*abar0[0]*favg[13]-0.2165063509461096*abar0[5]*favg[12]-0.2165063509461096*abar0[4]*favg[11]-0.375*abar0[2]*favg[7]-0.375*abar0[1]*favg[6]-0.2165063509461096*abar0[3]*favg[5]-0.375*abar0[0]*favg[3]-0.2165063509461096*abar0[2]*favg[2]-0.2165063509461096*abar0[1]*favg[1]-0.2165063509461096*abar0[0]*favg[0]; 
  incr[4] = 0.03227486121839514*B2[0]*favg[14]*dv2+0.08068715304598786*B2[0]*favg[13]*dv2+0.03608439182435161*B2[5]*favg[12]*dv2+0.03608439182435161*B2[4]*favg[11]*dv2+0.0625*B2[2]*favg[7]*dv2+0.0625*B2[1]*favg[6]*dv2+0.03608439182435161*B2[3]*favg[5]*dv2+0.0625*B2[0]*favg[3]*dv2+0.03608439182435161*favg[2]*B2[2]*dv2+0.03608439182435161*favg[1]*B2[1]*dv2+0.03608439182435161*favg[0]*B2[0]*dv2-0.4330127018922193*fjump[10]*amax-0.25*fjump[4]*amax+0.2165063509461096*abar0[0]*favg[10]+0.125*abar0[2]*favg[9]+0.125*abar0[1]*favg[8]+0.125*abar0[0]*favg[4]; 
  incr[5] = 0.0625*B2[3]*favg[10]*dv2+0.03608439182435161*B2[1]*favg[9]*dv2+0.03608439182435161*B2[2]*favg[8]*dv2+0.03608439182435161*B2[3]*favg[4]*dv2-0.25*fjump[5]*amax+0.2795084971874738*abar0[3]*favg[13]+0.1118033988749895*abar0[3]*favg[12]+0.1118033988749895*abar0[3]*favg[11]+0.2165063509461096*abar0[1]*favg[7]+0.2165063509461096*abar0[2]*favg[6]+0.1118033988749895*abar0[5]*favg[5]+0.1118033988749895*abar0[4]*favg[5]+0.125*abar0[0]*favg[5]+0.2165063509461096*abar0[3]*favg[3]+0.125*favg[0]*abar0[3]+0.125*abar0[1]*favg[2]+0.125*favg[1]*abar0[2]; 
  incr[6] = (-0.1082531754730548*B2[1]*favg[10]*dv2)-0.0625*B2[3]*favg[9]*dv2-0.05590169943749475*B2[4]*favg[8]*dv2-0.0625*B2[0]*favg[8]*dv2-0.0625*B2[1]*favg[4]*dv2+0.75*fjump[6]*amax+0.4330127018922193*fjump[1]*amax-0.4841229182759272*abar0[1]*favg[13]-0.1936491673103709*abar0[1]*favg[11]-0.375*abar0[3]*favg[7]-0.3354101966249685*abar0[4]*favg[6]-0.375*abar0[0]*favg[6]-0.2165063509461096*abar0[2]*favg[5]-0.1936491673103709*favg[1]*abar0[4]-0.375*abar0[1]*favg[3]-0.2165063509461096*favg[2]*abar0[3]-0.2165063509461096*abar0[0]*favg[1]-0.2165063509461096*favg[0]*abar0[1]; 
  incr[7] = (-0.1082531754730548*B2[2]*favg[10]*dv2)-0.05590169943749475*B2[5]*favg[9]*dv2-0.0625*B2[0]*favg[9]*dv2-0.0625*B2[3]*favg[8]*dv2-0.0625*B2[2]*favg[4]*dv2+0.75*fjump[7]*amax+0.4330127018922193*fjump[2]*amax-0.4841229182759272*abar0[2]*favg[13]-0.1936491673103709*abar0[2]*favg[12]-0.3354101966249685*abar0[5]*favg[7]-0.375*abar0[0]*favg[7]-0.375*abar0[3]*favg[6]-0.2165063509461096*abar0[1]*favg[5]-0.1936491673103709*favg[2]*abar0[5]-0.375*abar0[2]*favg[3]-0.2165063509461096*favg[1]*abar0[3]-0.2165063509461096*abar0[0]*favg[2]-0.2165063509461096*favg[0]*abar0[2]; 
  incr[8] = 0.03227486121839514*B2[1]*favg[14]*dv2+0.08068715304598786*B2[1]*favg[13]*dv2+0.03227486121839514*B2[1]*favg[11]*dv2+0.0625*B2[3]*favg[7]*dv2+0.05590169943749475*B2[4]*favg[6]*dv2+0.0625*B2[0]*favg[6]*dv2+0.03608439182435161*B2[2]*favg[5]*dv2+0.03227486121839514*favg[1]*B2[4]*dv2+0.03608439182435161*favg[2]*B2[3]*dv2+0.0625*B2[1]*favg[3]*dv2+0.03608439182435161*favg[0]*B2[1]*dv2+0.03608439182435161*B2[0]*favg[1]*dv2-0.25*fjump[8]*amax+0.2165063509461096*abar0[1]*favg[10]+0.125*abar0[3]*favg[9]+0.1118033988749895*abar0[4]*favg[8]+0.125*abar0[0]*favg[8]+0.125*abar0[1]*favg[4]; 
  incr[9] = 0.03227486121839514*B2[2]*favg[14]*dv2+0.08068715304598786*B2[2]*favg[13]*dv2+0.03227486121839514*B2[2]*favg[12]*dv2+0.05590169943749475*B2[5]*favg[7]*dv2+0.0625*B2[0]*favg[7]*dv2+0.0625*B2[3]*favg[6]*dv2+0.03227486121839514*favg[2]*B2[5]*dv2+0.03608439182435161*B2[1]*favg[5]*dv2+0.03608439182435161*favg[1]*B2[3]*dv2+0.0625*B2[2]*favg[3]*dv2+0.03608439182435161*favg[0]*B2[2]*dv2+0.03608439182435161*B2[0]*favg[2]*dv2-0.25*fjump[9]*amax+0.2165063509461096*abar0[2]*favg[10]+0.1118033988749895*abar0[5]*favg[9]+0.125*abar0[0]*favg[9]+0.125*abar0[3]*favg[8]+0.125*abar0[2]*favg[4]; 
  incr[10] = (-0.05590169943749475*B2[0]*favg[14]*dv2)-0.1397542485937369*B2[0]*favg[13]*dv2-0.0625*B2[5]*favg[12]*dv2-0.0625*B2[4]*favg[11]*dv2-0.1082531754730548*B2[2]*favg[7]*dv2-0.1082531754730548*B2[1]*favg[6]*dv2-0.0625*B2[3]*favg[5]*dv2-0.1082531754730548*B2[0]*favg[3]*dv2-0.0625*favg[2]*B2[2]*dv2-0.0625*favg[1]*B2[1]*dv2-0.0625*favg[0]*B2[0]*dv2+0.75*fjump[10]*amax+0.4330127018922193*fjump[4]*amax-0.375*abar0[0]*favg[10]-0.2165063509461096*abar0[2]*favg[9]-0.2165063509461096*abar0[1]*favg[8]-0.2165063509461096*abar0[0]*favg[4]; 
  incr[11] = 0.0625*B2[4]*favg[10]*dv2+0.03227486121839514*B2[1]*favg[8]*dv2+0.03608439182435162*favg[4]*B2[4]*dv2-0.25*fjump[11]*amax+0.2795084971874738*abar0[4]*favg[13]+0.0798595706249925*abar0[4]*favg[11]+0.125*abar0[0]*favg[11]+0.1936491673103709*abar0[1]*favg[6]+0.1118033988749895*abar0[3]*favg[5]+0.2165063509461096*favg[3]*abar0[4]+0.125*favg[0]*abar0[4]+0.1118033988749895*abar0[1]*favg[1]; 
  incr[12] = 0.0625*B2[5]*favg[10]*dv2+0.03227486121839514*B2[2]*favg[9]*dv2+0.03608439182435162*favg[4]*B2[5]*dv2-0.25*fjump[12]*amax+0.2795084971874738*abar0[5]*favg[13]+0.0798595706249925*abar0[5]*favg[12]+0.125*abar0[0]*favg[12]+0.1936491673103709*abar0[2]*favg[7]+0.1118033988749895*abar0[3]*favg[5]+0.2165063509461096*favg[3]*abar0[5]+0.125*favg[0]*abar0[5]+0.1118033988749895*abar0[2]*favg[2]; 
  incr[13] = 0.1397542485937369*B2[0]*favg[10]*dv2+0.08068715304598785*B2[2]*favg[9]*dv2+0.08068715304598785*B2[1]*favg[8]*dv2+0.08068715304598785*B2[0]*favg[4]*dv2-1.25*fjump[13]*amax-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[0]*amax+0.625*abar0[0]*favg[13]+0.2795084971874737*abar0[5]*favg[12]+0.2795084971874737*abar0[4]*favg[11]+0.484122918275927*abar0[2]*favg[7]+0.484122918275927*abar0[1]*favg[6]+0.2795084971874737*abar0[3]*favg[5]+0.484122918275927*abar0[0]*favg[3]+0.2795084971874737*abar0[2]*favg[2]+0.2795084971874737*abar0[1]*favg[1]+0.2795084971874737*abar0[0]*favg[0]; 
  incr[14] = 0.05590169943749475*B2[0]*favg[10]*dv2+0.03227486121839514*B2[2]*favg[9]*dv2+0.03227486121839514*B2[1]*favg[8]*dv2+0.03227486121839514*B2[0]*favg[4]*dv2-0.25*fjump[14]*amax+0.125*abar0[0]*favg[14]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 
  outr[5] += incr[5]*dv10; 
  outr[6] += incr[6]*dv10; 
  outr[7] += incr[7]*dv10; 
  outr[8] += incr[8]*dv10; 
  outr[9] += incr[9]*dv10; 
  outr[10] += incr[10]*dv10; 
  outr[11] += incr[11]*dv10; 
  outr[12] += incr[12]*dv10; 
  outr[13] += incr[13]*dv10; 
  outr[14] += incr[14]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += -1.0*incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += incr[6]*dv10; 
  outl[7] += incr[7]*dv10; 
  outl[8] += -1.0*incr[8]*dv10; 
  outl[9] += -1.0*incr[9]*dv10; 
  outl[10] += incr[10]*dv10; 
  outl[11] += -1.0*incr[11]*dv10; 
  outl[12] += -1.0*incr[12]*dv10; 
  outl[13] += -1.0*incr[13]*dv10; 
  outl[14] += -1.0*incr[14]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x2vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 


  abar1[0] = E1[0]-wv1*B2[0]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar1[2] = E1[2]-wv1*B2[2]; 

  double incr[5]; 

  double favg[5]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  double fjump[5]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = -1*fr[4]-fl[4]; 
  const double amid = 0.5*abar1[0]; 
  incr[0] = (-0.03608439182435161*B2[0]*favg[3]*dv1)-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.2165063509461096*abar1[0]*favg[4]+0.125*abar1[2]*favg[2]+0.125*abar1[1]*favg[1]+0.125*abar1[0]*favg[0]; 
  incr[1] = (-0.03608439182435161*B2[1]*favg[3]*dv1)-0.25*fjump[1]*amax+0.2165063509461096*abar1[1]*favg[4]+0.125*abar1[0]*favg[1]+0.125*favg[0]*abar1[1]; 
  incr[2] = (-0.03608439182435161*B2[2]*favg[3]*dv1)-0.25*fjump[2]*amax+0.2165063509461096*abar1[2]*favg[4]+0.125*abar1[0]*favg[2]+0.125*favg[0]*abar1[2]; 
  incr[3] = (-0.0625*B2[0]*favg[4]*dv1)-0.03608439182435161*favg[2]*B2[2]*dv1-0.03608439182435161*favg[1]*B2[1]*dv1-0.03608439182435161*favg[0]*B2[0]*dv1-0.25*fjump[3]*amax+0.125*abar1[0]*favg[3]; 
  incr[4] = 0.0625*B2[0]*favg[3]*dv1+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.375*abar1[0]*favg[4]-0.2165063509461096*abar1[2]*favg[2]-0.2165063509461096*abar1[1]*favg[1]-0.2165063509461096*abar1[0]*favg[0]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 
  outr[4] += incr[4]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += -1.0*incr[3]*dv11; 
  outl[4] += incr[4]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x2vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[6]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 


  abar1[0] = E1[0]-wv1*B2[0]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar1[2] = E1[2]-wv1*B2[2]; 
  abar1[3] = E1[3]-wv1*B2[3]; 
  abar1[4] = E1[4]-wv1*B2[4]; 
  abar1[5] = E1[5]-wv1*B2[5]; 

  double incr[15]; 

  double favg[15]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  double fjump[15]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = -1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = 1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = -1*fr[8]-fl[8]; 
  fjump[9] = -1*fr[9]-fl[9]; 
  fjump[10] = -1*fr[10]-fl[10]; 
  fjump[11] = 1*fr[11]-fl[11]; 
  fjump[12] = 1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = 1*fr[14]-fl[14]; 
  const double amid = (-0.5590169943749475*abar1[5])-0.5590169943749475*abar1[4]+0.5*abar1[0]; 
  incr[0] = (-0.0625*B2[0]*favg[10]*dv1)-0.03608439182435161*B2[2]*favg[7]*dv1-0.03608439182435161*B2[1]*favg[6]*dv1-0.03608439182435161*B2[0]*favg[3]*dv1-0.5590169943749475*fjump[14]*amax-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.2795084971874737*abar1[0]*favg[14]+0.125*abar1[5]*favg[12]+0.125*abar1[4]*favg[11]+0.2165063509461096*abar1[2]*favg[9]+0.2165063509461096*abar1[1]*favg[8]+0.125*abar1[3]*favg[5]+0.2165063509461096*abar1[0]*favg[4]+0.125*abar1[2]*favg[2]+0.125*abar1[1]*favg[1]+0.125*abar1[0]*favg[0]; 
  incr[1] = (-0.0625*B2[1]*favg[10]*dv1)-0.03608439182435161*B2[3]*favg[7]*dv1-0.03227486121839514*B2[4]*favg[6]*dv1-0.03608439182435161*B2[0]*favg[6]*dv1-0.03608439182435161*B2[1]*favg[3]*dv1-0.4330127018922193*fjump[8]*amax-0.25*fjump[1]*amax+0.2795084971874738*abar1[1]*favg[14]+0.1118033988749895*abar1[1]*favg[11]+0.2165063509461096*abar1[3]*favg[9]+0.1936491673103708*abar1[4]*favg[8]+0.2165063509461096*abar1[0]*favg[8]+0.125*abar1[2]*favg[5]+0.2165063509461096*abar1[1]*favg[4]+0.1118033988749895*favg[1]*abar1[4]+0.125*favg[2]*abar1[3]+0.125*abar1[0]*favg[1]+0.125*favg[0]*abar1[1]; 
  incr[2] = (-0.0625*B2[2]*favg[10]*dv1)-0.03227486121839514*B2[5]*favg[7]*dv1-0.03608439182435161*B2[0]*favg[7]*dv1-0.03608439182435161*B2[3]*favg[6]*dv1-0.03608439182435161*B2[2]*favg[3]*dv1-0.4330127018922193*fjump[9]*amax-0.25*fjump[2]*amax+0.2795084971874738*abar1[2]*favg[14]+0.1118033988749895*abar1[2]*favg[12]+0.1936491673103708*abar1[5]*favg[9]+0.2165063509461096*abar1[0]*favg[9]+0.2165063509461096*abar1[3]*favg[8]+0.125*abar1[1]*favg[5]+0.1118033988749895*favg[2]*abar1[5]+0.2165063509461096*abar1[2]*favg[4]+0.125*favg[1]*abar1[3]+0.125*abar1[0]*favg[2]+0.125*favg[0]*abar1[2]; 
  incr[3] = (-0.08068715304598786*B2[0]*favg[14]*dv1)-0.03227486121839514*B2[0]*favg[13]*dv1-0.03608439182435161*B2[5]*favg[12]*dv1-0.03608439182435161*B2[4]*favg[11]*dv1-0.0625*B2[2]*favg[9]*dv1-0.0625*B2[1]*favg[8]*dv1-0.03608439182435161*B2[3]*favg[5]*dv1-0.0625*B2[0]*favg[4]*dv1-0.03608439182435161*favg[2]*B2[2]*dv1-0.03608439182435161*favg[1]*B2[1]*dv1-0.03608439182435161*favg[0]*B2[0]*dv1-0.4330127018922193*fjump[10]*amax-0.25*fjump[3]*amax+0.2165063509461096*abar1[0]*favg[10]+0.125*abar1[2]*favg[7]+0.125*abar1[1]*favg[6]+0.125*abar1[0]*favg[3]; 
  incr[4] = 0.1082531754730548*B2[0]*favg[10]*dv1+0.0625*B2[2]*favg[7]*dv1+0.0625*B2[1]*favg[6]*dv1+0.0625*B2[0]*favg[3]*dv1+0.9682458365518543*fjump[14]*amax+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.4841229182759271*abar1[0]*favg[14]-0.2165063509461096*abar1[5]*favg[12]-0.2165063509461096*abar1[4]*favg[11]-0.375*abar1[2]*favg[9]-0.375*abar1[1]*favg[8]-0.2165063509461096*abar1[3]*favg[5]-0.375*abar1[0]*favg[4]-0.2165063509461096*abar1[2]*favg[2]-0.2165063509461096*abar1[1]*favg[1]-0.2165063509461096*abar1[0]*favg[0]; 
  incr[5] = (-0.0625*B2[3]*favg[10]*dv1)-0.03608439182435161*B2[1]*favg[7]*dv1-0.03608439182435161*B2[2]*favg[6]*dv1-0.03608439182435161*favg[3]*B2[3]*dv1-0.25*fjump[5]*amax+0.2795084971874738*abar1[3]*favg[14]+0.1118033988749895*abar1[3]*favg[12]+0.1118033988749895*abar1[3]*favg[11]+0.2165063509461096*abar1[1]*favg[9]+0.2165063509461096*abar1[2]*favg[8]+0.1118033988749895*abar1[5]*favg[5]+0.1118033988749895*abar1[4]*favg[5]+0.125*abar1[0]*favg[5]+0.2165063509461096*abar1[3]*favg[4]+0.125*favg[0]*abar1[3]+0.125*abar1[1]*favg[2]+0.125*favg[1]*abar1[2]; 
  incr[6] = (-0.08068715304598786*B2[1]*favg[14]*dv1)-0.03227486121839514*B2[1]*favg[13]*dv1-0.03227486121839514*B2[1]*favg[11]*dv1-0.0625*B2[3]*favg[9]*dv1-0.05590169943749475*B2[4]*favg[8]*dv1-0.0625*B2[0]*favg[8]*dv1-0.03608439182435161*B2[2]*favg[5]*dv1-0.03227486121839514*favg[1]*B2[4]*dv1-0.0625*B2[1]*favg[4]*dv1-0.03608439182435161*favg[2]*B2[3]*dv1-0.03608439182435161*favg[0]*B2[1]*dv1-0.03608439182435161*B2[0]*favg[1]*dv1-0.25*fjump[6]*amax+0.2165063509461096*abar1[1]*favg[10]+0.125*abar1[3]*favg[7]+0.1118033988749895*abar1[4]*favg[6]+0.125*abar1[0]*favg[6]+0.125*abar1[1]*favg[3]; 
  incr[7] = (-0.08068715304598786*B2[2]*favg[14]*dv1)-0.03227486121839514*B2[2]*favg[13]*dv1-0.03227486121839514*B2[2]*favg[12]*dv1-0.05590169943749475*B2[5]*favg[9]*dv1-0.0625*B2[0]*favg[9]*dv1-0.0625*B2[3]*favg[8]*dv1-0.03227486121839514*favg[2]*B2[5]*dv1-0.03608439182435161*B2[1]*favg[5]*dv1-0.0625*B2[2]*favg[4]*dv1-0.03608439182435161*favg[1]*B2[3]*dv1-0.03608439182435161*favg[0]*B2[2]*dv1-0.03608439182435161*B2[0]*favg[2]*dv1-0.25*fjump[7]*amax+0.2165063509461096*abar1[2]*favg[10]+0.1118033988749895*abar1[5]*favg[7]+0.125*abar1[0]*favg[7]+0.125*abar1[3]*favg[6]+0.125*abar1[2]*favg[3]; 
  incr[8] = 0.1082531754730548*B2[1]*favg[10]*dv1+0.0625*B2[3]*favg[7]*dv1+0.05590169943749475*B2[4]*favg[6]*dv1+0.0625*B2[0]*favg[6]*dv1+0.0625*B2[1]*favg[3]*dv1+0.75*fjump[8]*amax+0.4330127018922193*fjump[1]*amax-0.4841229182759272*abar1[1]*favg[14]-0.1936491673103709*abar1[1]*favg[11]-0.375*abar1[3]*favg[9]-0.3354101966249685*abar1[4]*favg[8]-0.375*abar1[0]*favg[8]-0.2165063509461096*abar1[2]*favg[5]-0.375*abar1[1]*favg[4]-0.1936491673103709*favg[1]*abar1[4]-0.2165063509461096*favg[2]*abar1[3]-0.2165063509461096*abar1[0]*favg[1]-0.2165063509461096*favg[0]*abar1[1]; 
  incr[9] = 0.1082531754730548*B2[2]*favg[10]*dv1+0.05590169943749475*B2[5]*favg[7]*dv1+0.0625*B2[0]*favg[7]*dv1+0.0625*B2[3]*favg[6]*dv1+0.0625*B2[2]*favg[3]*dv1+0.75*fjump[9]*amax+0.4330127018922193*fjump[2]*amax-0.4841229182759272*abar1[2]*favg[14]-0.1936491673103709*abar1[2]*favg[12]-0.3354101966249685*abar1[5]*favg[9]-0.375*abar1[0]*favg[9]-0.375*abar1[3]*favg[8]-0.2165063509461096*abar1[1]*favg[5]-0.1936491673103709*favg[2]*abar1[5]-0.375*abar1[2]*favg[4]-0.2165063509461096*favg[1]*abar1[3]-0.2165063509461096*abar1[0]*favg[2]-0.2165063509461096*favg[0]*abar1[2]; 
  incr[10] = 0.1397542485937369*B2[0]*favg[14]*dv1+0.05590169943749475*B2[0]*favg[13]*dv1+0.0625*B2[5]*favg[12]*dv1+0.0625*B2[4]*favg[11]*dv1+0.1082531754730548*B2[2]*favg[9]*dv1+0.1082531754730548*B2[1]*favg[8]*dv1+0.0625*B2[3]*favg[5]*dv1+0.1082531754730548*B2[0]*favg[4]*dv1+0.0625*favg[2]*B2[2]*dv1+0.0625*favg[1]*B2[1]*dv1+0.0625*favg[0]*B2[0]*dv1+0.75*fjump[10]*amax+0.4330127018922193*fjump[3]*amax-0.375*abar1[0]*favg[10]-0.2165063509461096*abar1[2]*favg[7]-0.2165063509461096*abar1[1]*favg[6]-0.2165063509461096*abar1[0]*favg[3]; 
  incr[11] = (-0.0625*B2[4]*favg[10]*dv1)-0.03227486121839514*B2[1]*favg[6]*dv1-0.03608439182435162*favg[3]*B2[4]*dv1-0.25*fjump[11]*amax+0.2795084971874738*abar1[4]*favg[14]+0.0798595706249925*abar1[4]*favg[11]+0.125*abar1[0]*favg[11]+0.1936491673103709*abar1[1]*favg[8]+0.1118033988749895*abar1[3]*favg[5]+0.2165063509461096*abar1[4]*favg[4]+0.125*favg[0]*abar1[4]+0.1118033988749895*abar1[1]*favg[1]; 
  incr[12] = (-0.0625*B2[5]*favg[10]*dv1)-0.03227486121839514*B2[2]*favg[7]*dv1-0.03608439182435162*favg[3]*B2[5]*dv1-0.25*fjump[12]*amax+0.2795084971874738*abar1[5]*favg[14]+0.0798595706249925*abar1[5]*favg[12]+0.125*abar1[0]*favg[12]+0.1936491673103709*abar1[2]*favg[9]+0.1118033988749895*abar1[3]*favg[5]+0.2165063509461096*favg[4]*abar1[5]+0.125*favg[0]*abar1[5]+0.1118033988749895*abar1[2]*favg[2]; 
  incr[13] = (-0.05590169943749475*B2[0]*favg[10]*dv1)-0.03227486121839514*B2[2]*favg[7]*dv1-0.03227486121839514*B2[1]*favg[6]*dv1-0.03227486121839514*B2[0]*favg[3]*dv1-0.25*fjump[13]*amax+0.125*abar1[0]*favg[13]; 
  incr[14] = (-0.1397542485937369*B2[0]*favg[10]*dv1)-0.08068715304598785*B2[2]*favg[7]*dv1-0.08068715304598785*B2[1]*favg[6]*dv1-0.08068715304598785*B2[0]*favg[3]*dv1-1.25*fjump[14]*amax-0.9682458365518543*fjump[4]*amax-0.5590169943749475*fjump[0]*amax+0.625*abar1[0]*favg[14]+0.2795084971874737*abar1[5]*favg[12]+0.2795084971874737*abar1[4]*favg[11]+0.484122918275927*abar1[2]*favg[9]+0.484122918275927*abar1[1]*favg[8]+0.2795084971874737*abar1[3]*favg[5]+0.484122918275927*abar1[0]*favg[4]+0.2795084971874737*abar1[2]*favg[2]+0.2795084971874737*abar1[1]*favg[1]+0.2795084971874737*abar1[0]*favg[0]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 
  outr[4] += incr[4]*dv11; 
  outr[5] += incr[5]*dv11; 
  outr[6] += incr[6]*dv11; 
  outr[7] += incr[7]*dv11; 
  outr[8] += incr[8]*dv11; 
  outr[9] += incr[9]*dv11; 
  outr[10] += incr[10]*dv11; 
  outr[11] += incr[11]*dv11; 
  outr[12] += incr[12]*dv11; 
  outr[13] += incr[13]*dv11; 
  outr[14] += incr[14]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += -1.0*incr[3]*dv11; 
  outl[4] += incr[4]*dv11; 
  outl[5] += -1.0*incr[5]*dv11; 
  outl[6] += -1.0*incr[6]*dv11; 
  outl[7] += -1.0*incr[7]*dv11; 
  outl[8] += incr[8]*dv11; 
  outl[9] += incr[9]*dv11; 
  outl[10] += incr[10]*dv11; 
  outl[11] += -1.0*incr[11]*dv11; 
  outl[12] += -1.0*incr[12]*dv11; 
  outl[13] += -1.0*incr[13]*dv11; 
  outl[14] += -1.0*incr[14]*dv11; 
return std::abs(amid); 
} 
