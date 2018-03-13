#include <VlasovModDecl.h> 
double VlasovSurfElcMag2x3vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar0[0] = E0[0]+wv2*B2[0]-wv3*B1[0]; 
  abar0[1] = E0[1]+wv2*B2[1]-wv3*B1[1]; 
  abar0[2] = E0[2]+wv2*B2[2]-wv3*B1[2]; 

  double incr[6]; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  const double amid = 0.5*abar0[0]; 
  incr[0] = (-0.03608439182435161*B1[0]*favg[5]*dv3)+0.03608439182435161*B2[0]*favg[4]*dv2-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.2165063509461096*abar0[0]*favg[3]+0.125*abar0[2]*favg[2]+0.125*abar0[1]*favg[1]+0.125*abar0[0]*favg[0]; 
  incr[1] = (-0.03608439182435161*B1[1]*favg[5]*dv3)+0.03608439182435161*B2[1]*favg[4]*dv2-0.25*fjump[1]*amax+0.2165063509461096*abar0[1]*favg[3]+0.125*abar0[0]*favg[1]+0.125*favg[0]*abar0[1]; 
  incr[2] = (-0.03608439182435161*B1[2]*favg[5]*dv3)+0.03608439182435161*B2[2]*favg[4]*dv2-0.25*fjump[2]*amax+0.2165063509461096*abar0[2]*favg[3]+0.125*abar0[0]*favg[2]+0.125*favg[0]*abar0[2]; 
  incr[3] = 0.0625*B1[0]*favg[5]*dv3-0.0625*B2[0]*favg[4]*dv2+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.375*abar0[0]*favg[3]-0.2165063509461096*abar0[2]*favg[2]-0.2165063509461096*abar0[1]*favg[1]-0.2165063509461096*abar0[0]*favg[0]; 
  incr[4] = 0.0625*B2[0]*favg[3]*dv2+0.03608439182435161*favg[2]*B2[2]*dv2+0.03608439182435161*favg[1]*B2[1]*dv2+0.03608439182435161*favg[0]*B2[0]*dv2-0.25*fjump[4]*amax+0.125*abar0[0]*favg[4]; 
  incr[5] = (-0.0625*B1[0]*favg[3]*dv3)-0.03608439182435161*favg[2]*B1[2]*dv3-0.03608439182435161*favg[1]*B1[1]*dv3-0.03608439182435161*favg[0]*B1[0]*dv3-0.25*fjump[5]*amax+0.125*abar0[0]*favg[5]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 
  outr[5] += incr[5]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += -1.0*incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 

  double abar2[6]; 


  abar0[0] = E0[0]+wv2*B2[0]-wv3*B1[0]; 
  abar0[1] = E0[1]+wv2*B2[1]-wv3*B1[1]; 
  abar0[2] = E0[2]+wv2*B2[2]-wv3*B1[2]; 
  abar0[3] = E0[3]+wv2*B2[3]-wv3*B1[3]; 
  abar0[4] = E0[4]+wv2*B2[4]-wv3*B1[4]; 
  abar0[5] = E0[5]+wv2*B2[5]-wv3*B1[5]; 

  double incr[21]; 

  double favg[21]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = -1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = -1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = 1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = -1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = 1*fr[6]-fl[6]; 
  fjump[7] = -1*fr[7]-fl[7]; 
  fjump[8] = -1*fr[8]-fl[8]; 
  fjump[9] = 1*fr[9]-fl[9]; 
  fjump[10] = 1*fr[10]-fl[10]; 
  fjump[11] = -1*fr[11]-fl[11]; 
  fjump[12] = 1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = -1*fr[14]-fl[14]; 
  fjump[15] = 1*fr[15]-fl[15]; 
  fjump[16] = 1*fr[16]-fl[16]; 
  fjump[17] = 1*fr[17]-fl[17]; 
  fjump[18] = 1*fr[18]-fl[18]; 
  fjump[19] = 1*fr[19]-fl[19]; 
  fjump[20] = 1*fr[20]-fl[20]; 
  const double amid = (-0.5590169943749475*abar0[5])-0.5590169943749475*abar0[4]+0.5*abar0[0]; 
  incr[0] = (-0.0625*B1[0]*favg[14]*dv3)-0.03608439182435161*B1[2]*favg[13]*dv3-0.03608439182435161*B1[1]*favg[12]*dv3-0.03608439182435161*B1[0]*favg[5]*dv3+0.0625*B2[0]*favg[11]*dv2+0.03608439182435161*B2[2]*favg[10]*dv2+0.03608439182435161*B2[1]*favg[9]*dv2+0.03608439182435161*B2[0]*favg[4]*dv2-0.5590169943749475*fjump[18]*amax-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.2795084971874737*abar0[0]*favg[18]+0.125*abar0[5]*favg[17]+0.125*abar0[4]*favg[16]+0.2165063509461096*abar0[2]*favg[8]+0.2165063509461096*abar0[1]*favg[7]+0.125*abar0[3]*favg[6]+0.2165063509461096*abar0[0]*favg[3]+0.125*abar0[2]*favg[2]+0.125*abar0[1]*favg[1]+0.125*abar0[0]*favg[0]; 
  incr[1] = (-0.0625*B1[1]*favg[14]*dv3)-0.03608439182435161*B1[3]*favg[13]*dv3-0.03227486121839514*B1[4]*favg[12]*dv3-0.03608439182435161*B1[0]*favg[12]*dv3-0.03608439182435161*B1[1]*favg[5]*dv3+0.0625*B2[1]*favg[11]*dv2+0.03608439182435161*B2[3]*favg[10]*dv2+0.03227486121839514*B2[4]*favg[9]*dv2+0.03608439182435161*B2[0]*favg[9]*dv2+0.03608439182435161*B2[1]*favg[4]*dv2-0.4330127018922193*fjump[7]*amax-0.25*fjump[1]*amax+0.2795084971874738*abar0[1]*favg[18]+0.1118033988749895*abar0[1]*favg[16]+0.2165063509461096*abar0[3]*favg[8]+0.1936491673103708*abar0[4]*favg[7]+0.2165063509461096*abar0[0]*favg[7]+0.125*abar0[2]*favg[6]+0.1118033988749895*favg[1]*abar0[4]+0.2165063509461096*abar0[1]*favg[3]+0.125*favg[2]*abar0[3]+0.125*abar0[0]*favg[1]+0.125*favg[0]*abar0[1]; 
  incr[2] = (-0.0625*B1[2]*favg[14]*dv3)-0.03227486121839514*B1[5]*favg[13]*dv3-0.03608439182435161*B1[0]*favg[13]*dv3-0.03608439182435161*B1[3]*favg[12]*dv3-0.03608439182435161*B1[2]*favg[5]*dv3+0.0625*B2[2]*favg[11]*dv2+0.03227486121839514*B2[5]*favg[10]*dv2+0.03608439182435161*B2[0]*favg[10]*dv2+0.03608439182435161*B2[3]*favg[9]*dv2+0.03608439182435161*B2[2]*favg[4]*dv2-0.4330127018922193*fjump[8]*amax-0.25*fjump[2]*amax+0.2795084971874738*abar0[2]*favg[18]+0.1118033988749895*abar0[2]*favg[17]+0.1936491673103708*abar0[5]*favg[8]+0.2165063509461096*abar0[0]*favg[8]+0.2165063509461096*abar0[3]*favg[7]+0.125*abar0[1]*favg[6]+0.1118033988749895*favg[2]*abar0[5]+0.2165063509461096*abar0[2]*favg[3]+0.125*favg[1]*abar0[3]+0.125*abar0[0]*favg[2]+0.125*favg[0]*abar0[2]; 
  incr[3] = 0.1082531754730548*B1[0]*favg[14]*dv3+0.0625*B1[2]*favg[13]*dv3+0.0625*B1[1]*favg[12]*dv3+0.0625*B1[0]*favg[5]*dv3-0.1082531754730548*B2[0]*favg[11]*dv2-0.0625*B2[2]*favg[10]*dv2-0.0625*B2[1]*favg[9]*dv2-0.0625*B2[0]*favg[4]*dv2+0.9682458365518543*fjump[18]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.4841229182759271*abar0[0]*favg[18]-0.2165063509461096*abar0[5]*favg[17]-0.2165063509461096*abar0[4]*favg[16]-0.375*abar0[2]*favg[8]-0.375*abar0[1]*favg[7]-0.2165063509461096*abar0[3]*favg[6]-0.375*abar0[0]*favg[3]-0.2165063509461096*abar0[2]*favg[2]-0.2165063509461096*abar0[1]*favg[1]-0.2165063509461096*abar0[0]*favg[0]; 
  incr[4] = (-0.03608439182435161*B1[0]*favg[15]*dv3)+0.03227486121839514*B2[0]*favg[19]*dv2+0.08068715304598786*B2[0]*favg[18]*dv2+0.03608439182435161*B2[5]*favg[17]*dv2+0.03608439182435161*B2[4]*favg[16]*dv2+0.0625*B2[2]*favg[8]*dv2+0.0625*B2[1]*favg[7]*dv2+0.03608439182435161*B2[3]*favg[6]*dv2+0.0625*B2[0]*favg[3]*dv2+0.03608439182435161*favg[2]*B2[2]*dv2+0.03608439182435161*favg[1]*B2[1]*dv2+0.03608439182435161*favg[0]*B2[0]*dv2-0.4330127018922193*fjump[11]*amax-0.25*fjump[4]*amax+0.2165063509461096*abar0[0]*favg[11]+0.125*abar0[2]*favg[10]+0.125*abar0[1]*favg[9]+0.125*abar0[0]*favg[4]; 
  incr[5] = (-0.03227486121839514*B1[0]*favg[20]*dv3)-0.08068715304598786*B1[0]*favg[18]*dv3-0.03608439182435161*B1[5]*favg[17]*dv3-0.03608439182435161*B1[4]*favg[16]*dv3-0.0625*B1[2]*favg[8]*dv3-0.0625*B1[1]*favg[7]*dv3-0.03608439182435161*B1[3]*favg[6]*dv3-0.0625*B1[0]*favg[3]*dv3-0.03608439182435161*favg[2]*B1[2]*dv3-0.03608439182435161*favg[1]*B1[1]*dv3-0.03608439182435161*favg[0]*B1[0]*dv3+0.03608439182435161*B2[0]*favg[15]*dv2-0.4330127018922193*fjump[14]*amax-0.25*fjump[5]*amax+0.2165063509461096*abar0[0]*favg[14]+0.125*abar0[2]*favg[13]+0.125*abar0[1]*favg[12]+0.125*abar0[0]*favg[5]; 
  incr[6] = (-0.0625*B1[3]*favg[14]*dv3)-0.03608439182435161*B1[1]*favg[13]*dv3-0.03608439182435161*B1[2]*favg[12]*dv3-0.03608439182435161*B1[3]*favg[5]*dv3+0.0625*B2[3]*favg[11]*dv2+0.03608439182435161*B2[1]*favg[10]*dv2+0.03608439182435161*B2[2]*favg[9]*dv2+0.03608439182435161*B2[3]*favg[4]*dv2-0.25*fjump[6]*amax+0.2795084971874738*abar0[3]*favg[18]+0.1118033988749895*abar0[3]*favg[17]+0.1118033988749895*abar0[3]*favg[16]+0.2165063509461096*abar0[1]*favg[8]+0.2165063509461096*abar0[2]*favg[7]+0.1118033988749895*abar0[5]*favg[6]+0.1118033988749895*abar0[4]*favg[6]+0.125*abar0[0]*favg[6]+0.2165063509461096*abar0[3]*favg[3]+0.125*favg[0]*abar0[3]+0.125*abar0[1]*favg[2]+0.125*favg[1]*abar0[2]; 
  incr[7] = 0.1082531754730548*B1[1]*favg[14]*dv3+0.0625*B1[3]*favg[13]*dv3+0.05590169943749475*B1[4]*favg[12]*dv3+0.0625*B1[0]*favg[12]*dv3+0.0625*B1[1]*favg[5]*dv3-0.1082531754730548*B2[1]*favg[11]*dv2-0.0625*B2[3]*favg[10]*dv2-0.05590169943749475*B2[4]*favg[9]*dv2-0.0625*B2[0]*favg[9]*dv2-0.0625*B2[1]*favg[4]*dv2+0.75*fjump[7]*amax+0.4330127018922193*fjump[1]*amax-0.4841229182759272*abar0[1]*favg[18]-0.1936491673103709*abar0[1]*favg[16]-0.375*abar0[3]*favg[8]-0.3354101966249685*abar0[4]*favg[7]-0.375*abar0[0]*favg[7]-0.2165063509461096*abar0[2]*favg[6]-0.1936491673103709*favg[1]*abar0[4]-0.375*abar0[1]*favg[3]-0.2165063509461096*favg[2]*abar0[3]-0.2165063509461096*abar0[0]*favg[1]-0.2165063509461096*favg[0]*abar0[1]; 
  incr[8] = 0.1082531754730548*B1[2]*favg[14]*dv3+0.05590169943749475*B1[5]*favg[13]*dv3+0.0625*B1[0]*favg[13]*dv3+0.0625*B1[3]*favg[12]*dv3+0.0625*B1[2]*favg[5]*dv3-0.1082531754730548*B2[2]*favg[11]*dv2-0.05590169943749475*B2[5]*favg[10]*dv2-0.0625*B2[0]*favg[10]*dv2-0.0625*B2[3]*favg[9]*dv2-0.0625*B2[2]*favg[4]*dv2+0.75*fjump[8]*amax+0.4330127018922193*fjump[2]*amax-0.4841229182759272*abar0[2]*favg[18]-0.1936491673103709*abar0[2]*favg[17]-0.3354101966249685*abar0[5]*favg[8]-0.375*abar0[0]*favg[8]-0.375*abar0[3]*favg[7]-0.2165063509461096*abar0[1]*favg[6]-0.1936491673103709*favg[2]*abar0[5]-0.375*abar0[2]*favg[3]-0.2165063509461096*favg[1]*abar0[3]-0.2165063509461096*abar0[0]*favg[2]-0.2165063509461096*favg[0]*abar0[2]; 
  incr[9] = (-0.03608439182435161*B1[1]*favg[15]*dv3)+0.03227486121839514*B2[1]*favg[19]*dv2+0.08068715304598786*B2[1]*favg[18]*dv2+0.03227486121839514*B2[1]*favg[16]*dv2+0.0625*B2[3]*favg[8]*dv2+0.05590169943749475*B2[4]*favg[7]*dv2+0.0625*B2[0]*favg[7]*dv2+0.03608439182435161*B2[2]*favg[6]*dv2+0.03227486121839514*favg[1]*B2[4]*dv2+0.03608439182435161*favg[2]*B2[3]*dv2+0.0625*B2[1]*favg[3]*dv2+0.03608439182435161*favg[0]*B2[1]*dv2+0.03608439182435161*B2[0]*favg[1]*dv2-0.25*fjump[9]*amax+0.2165063509461096*abar0[1]*favg[11]+0.125*abar0[3]*favg[10]+0.1118033988749895*abar0[4]*favg[9]+0.125*abar0[0]*favg[9]+0.125*abar0[1]*favg[4]; 
  incr[10] = (-0.03608439182435161*B1[2]*favg[15]*dv3)+0.03227486121839514*B2[2]*favg[19]*dv2+0.08068715304598786*B2[2]*favg[18]*dv2+0.03227486121839514*B2[2]*favg[17]*dv2+0.05590169943749475*B2[5]*favg[8]*dv2+0.0625*B2[0]*favg[8]*dv2+0.0625*B2[3]*favg[7]*dv2+0.03608439182435161*B2[1]*favg[6]*dv2+0.03227486121839514*favg[2]*B2[5]*dv2+0.03608439182435161*favg[1]*B2[3]*dv2+0.0625*B2[2]*favg[3]*dv2+0.03608439182435161*favg[0]*B2[2]*dv2+0.03608439182435161*B2[0]*favg[2]*dv2-0.25*fjump[10]*amax+0.2165063509461096*abar0[2]*favg[11]+0.1118033988749895*abar0[5]*favg[10]+0.125*abar0[0]*favg[10]+0.125*abar0[3]*favg[9]+0.125*abar0[2]*favg[4]; 
  incr[11] = 0.0625*B1[0]*favg[15]*dv3-0.05590169943749475*B2[0]*favg[19]*dv2-0.1397542485937369*B2[0]*favg[18]*dv2-0.0625*B2[5]*favg[17]*dv2-0.0625*B2[4]*favg[16]*dv2-0.1082531754730548*B2[2]*favg[8]*dv2-0.1082531754730548*B2[1]*favg[7]*dv2-0.0625*B2[3]*favg[6]*dv2-0.1082531754730548*B2[0]*favg[3]*dv2-0.0625*favg[2]*B2[2]*dv2-0.0625*favg[1]*B2[1]*dv2-0.0625*favg[0]*B2[0]*dv2+0.75*fjump[11]*amax+0.4330127018922193*fjump[4]*amax-0.375*abar0[0]*favg[11]-0.2165063509461096*abar0[2]*favg[10]-0.2165063509461096*abar0[1]*favg[9]-0.2165063509461096*abar0[0]*favg[4]; 
  incr[12] = (-0.03227486121839514*B1[1]*favg[20]*dv3)-0.08068715304598786*B1[1]*favg[18]*dv3-0.03227486121839514*B1[1]*favg[16]*dv3-0.0625*B1[3]*favg[8]*dv3-0.05590169943749475*B1[4]*favg[7]*dv3-0.0625*B1[0]*favg[7]*dv3-0.03608439182435161*B1[2]*favg[6]*dv3-0.03227486121839514*favg[1]*B1[4]*dv3-0.03608439182435161*favg[2]*B1[3]*dv3-0.0625*B1[1]*favg[3]*dv3-0.03608439182435161*favg[0]*B1[1]*dv3-0.03608439182435161*B1[0]*favg[1]*dv3+0.03608439182435161*B2[1]*favg[15]*dv2-0.25*fjump[12]*amax+0.2165063509461096*abar0[1]*favg[14]+0.125*abar0[3]*favg[13]+0.1118033988749895*abar0[4]*favg[12]+0.125*abar0[0]*favg[12]+0.125*abar0[1]*favg[5]; 
  incr[13] = (-0.03227486121839514*B1[2]*favg[20]*dv3)-0.08068715304598786*B1[2]*favg[18]*dv3-0.03227486121839514*B1[2]*favg[17]*dv3-0.05590169943749475*B1[5]*favg[8]*dv3-0.0625*B1[0]*favg[8]*dv3-0.0625*B1[3]*favg[7]*dv3-0.03608439182435161*B1[1]*favg[6]*dv3-0.03227486121839514*favg[2]*B1[5]*dv3-0.03608439182435161*favg[1]*B1[3]*dv3-0.0625*B1[2]*favg[3]*dv3-0.03608439182435161*favg[0]*B1[2]*dv3-0.03608439182435161*B1[0]*favg[2]*dv3+0.03608439182435161*B2[2]*favg[15]*dv2-0.25*fjump[13]*amax+0.2165063509461096*abar0[2]*favg[14]+0.1118033988749895*abar0[5]*favg[13]+0.125*abar0[0]*favg[13]+0.125*abar0[3]*favg[12]+0.125*abar0[2]*favg[5]; 
  incr[14] = 0.05590169943749475*B1[0]*favg[20]*dv3+0.1397542485937369*B1[0]*favg[18]*dv3+0.0625*B1[5]*favg[17]*dv3+0.0625*B1[4]*favg[16]*dv3+0.1082531754730548*B1[2]*favg[8]*dv3+0.1082531754730548*B1[1]*favg[7]*dv3+0.0625*B1[3]*favg[6]*dv3+0.1082531754730548*B1[0]*favg[3]*dv3+0.0625*favg[2]*B1[2]*dv3+0.0625*favg[1]*B1[1]*dv3+0.0625*favg[0]*B1[0]*dv3-0.0625*B2[0]*favg[15]*dv2+0.75*fjump[14]*amax+0.4330127018922193*fjump[5]*amax-0.375*abar0[0]*favg[14]-0.2165063509461096*abar0[2]*favg[13]-0.2165063509461096*abar0[1]*favg[12]-0.2165063509461096*abar0[0]*favg[5]; 
  incr[15] = (-0.0625*B1[0]*favg[11]*dv3)-0.03608439182435161*B1[2]*favg[10]*dv3-0.03608439182435161*B1[1]*favg[9]*dv3-0.03608439182435161*B1[0]*favg[4]*dv3+0.0625*B2[0]*favg[14]*dv2+0.03608439182435161*B2[2]*favg[13]*dv2+0.03608439182435161*B2[1]*favg[12]*dv2+0.03608439182435161*B2[0]*favg[5]*dv2-0.25*fjump[15]*amax+0.125*abar0[0]*favg[15]; 
  incr[16] = (-0.0625*B1[4]*favg[14]*dv3)-0.03227486121839514*B1[1]*favg[12]*dv3-0.03608439182435162*B1[4]*favg[5]*dv3+0.0625*B2[4]*favg[11]*dv2+0.03227486121839514*B2[1]*favg[9]*dv2+0.03608439182435162*favg[4]*B2[4]*dv2-0.25*fjump[16]*amax+0.2795084971874738*abar0[4]*favg[18]+0.0798595706249925*abar0[4]*favg[16]+0.125*abar0[0]*favg[16]+0.1936491673103709*abar0[1]*favg[7]+0.1118033988749895*abar0[3]*favg[6]+0.2165063509461096*favg[3]*abar0[4]+0.125*favg[0]*abar0[4]+0.1118033988749895*abar0[1]*favg[1]; 
  incr[17] = (-0.0625*B1[5]*favg[14]*dv3)-0.03227486121839514*B1[2]*favg[13]*dv3-0.03608439182435162*favg[5]*B1[5]*dv3+0.0625*B2[5]*favg[11]*dv2+0.03227486121839514*B2[2]*favg[10]*dv2+0.03608439182435162*favg[4]*B2[5]*dv2-0.25*fjump[17]*amax+0.2795084971874738*abar0[5]*favg[18]+0.0798595706249925*abar0[5]*favg[17]+0.125*abar0[0]*favg[17]+0.1936491673103709*abar0[2]*favg[8]+0.1118033988749895*abar0[3]*favg[6]+0.2165063509461096*favg[3]*abar0[5]+0.125*favg[0]*abar0[5]+0.1118033988749895*abar0[2]*favg[2]; 
  incr[18] = (-0.1397542485937369*B1[0]*favg[14]*dv3)-0.08068715304598785*B1[2]*favg[13]*dv3-0.08068715304598785*B1[1]*favg[12]*dv3-0.08068715304598785*B1[0]*favg[5]*dv3+0.1397542485937369*B2[0]*favg[11]*dv2+0.08068715304598785*B2[2]*favg[10]*dv2+0.08068715304598785*B2[1]*favg[9]*dv2+0.08068715304598785*B2[0]*favg[4]*dv2-1.25*fjump[18]*amax-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[0]*amax+0.625*abar0[0]*favg[18]+0.2795084971874737*abar0[5]*favg[17]+0.2795084971874737*abar0[4]*favg[16]+0.484122918275927*abar0[2]*favg[8]+0.484122918275927*abar0[1]*favg[7]+0.2795084971874737*abar0[3]*favg[6]+0.484122918275927*abar0[0]*favg[3]+0.2795084971874737*abar0[2]*favg[2]+0.2795084971874737*abar0[1]*favg[1]+0.2795084971874737*abar0[0]*favg[0]; 
  incr[19] = 0.05590169943749475*B2[0]*favg[11]*dv2+0.03227486121839514*B2[2]*favg[10]*dv2+0.03227486121839514*B2[1]*favg[9]*dv2+0.03227486121839514*B2[0]*favg[4]*dv2-0.25*fjump[19]*amax+0.125*abar0[0]*favg[19]; 
  incr[20] = (-0.05590169943749475*B1[0]*favg[14]*dv3)-0.03227486121839514*B1[2]*favg[13]*dv3-0.03227486121839514*B1[1]*favg[12]*dv3-0.03227486121839514*B1[0]*favg[5]*dv3-0.25*fjump[20]*amax+0.125*abar0[0]*favg[20]; 

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
  outr[15] += incr[15]*dv10; 
  outr[16] += incr[16]*dv10; 
  outr[17] += incr[17]*dv10; 
  outr[18] += incr[18]*dv10; 
  outr[19] += incr[19]*dv10; 
  outr[20] += incr[20]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += -1.0*incr[2]*dv10; 
  outl[3] += incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += -1.0*incr[5]*dv10; 
  outl[6] += -1.0*incr[6]*dv10; 
  outl[7] += incr[7]*dv10; 
  outl[8] += incr[8]*dv10; 
  outl[9] += -1.0*incr[9]*dv10; 
  outl[10] += -1.0*incr[10]*dv10; 
  outl[11] += incr[11]*dv10; 
  outl[12] += -1.0*incr[12]*dv10; 
  outl[13] += -1.0*incr[13]*dv10; 
  outl[14] += incr[14]*dv10; 
  outl[15] += -1.0*incr[15]*dv10; 
  outl[16] += -1.0*incr[16]*dv10; 
  outl[17] += -1.0*incr[17]*dv10; 
  outl[18] += -1.0*incr[18]*dv10; 
  outl[19] += -1.0*incr[19]*dv10; 
  outl[20] += -1.0*incr[20]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar1[0] = E1[0]+wv3*B0[0]-wv1*B2[0]; 
  abar1[1] = E1[1]+wv3*B0[1]-wv1*B2[1]; 
  abar1[2] = E1[2]+wv3*B0[2]-wv1*B2[2]; 

  double incr[6]; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = -1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  const double amid = 0.5*abar1[0]; 
  incr[0] = 0.03608439182435161*B0[0]*favg[5]*dv3-0.03608439182435161*B2[0]*favg[3]*dv1-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.2165063509461096*abar1[0]*favg[4]+0.125*abar1[2]*favg[2]+0.125*abar1[1]*favg[1]+0.125*abar1[0]*favg[0]; 
  incr[1] = 0.03608439182435161*B0[1]*favg[5]*dv3-0.03608439182435161*B2[1]*favg[3]*dv1-0.25*fjump[1]*amax+0.2165063509461096*abar1[1]*favg[4]+0.125*abar1[0]*favg[1]+0.125*favg[0]*abar1[1]; 
  incr[2] = 0.03608439182435161*B0[2]*favg[5]*dv3-0.03608439182435161*B2[2]*favg[3]*dv1-0.25*fjump[2]*amax+0.2165063509461096*abar1[2]*favg[4]+0.125*abar1[0]*favg[2]+0.125*favg[0]*abar1[2]; 
  incr[3] = (-0.0625*B2[0]*favg[4]*dv1)-0.03608439182435161*favg[2]*B2[2]*dv1-0.03608439182435161*favg[1]*B2[1]*dv1-0.03608439182435161*favg[0]*B2[0]*dv1-0.25*fjump[3]*amax+0.125*abar1[0]*favg[3]; 
  incr[4] = (-0.0625*B0[0]*favg[5]*dv3)+0.0625*B2[0]*favg[3]*dv1+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.375*abar1[0]*favg[4]-0.2165063509461096*abar1[2]*favg[2]-0.2165063509461096*abar1[1]*favg[1]-0.2165063509461096*abar1[0]*favg[0]; 
  incr[5] = 0.0625*B0[0]*favg[4]*dv3+0.03608439182435161*favg[2]*B0[2]*dv3+0.03608439182435161*favg[1]*B0[1]*dv3+0.03608439182435161*favg[0]*B0[0]*dv3-0.25*fjump[5]*amax+0.125*abar1[0]*favg[5]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 
  outr[4] += incr[4]*dv11; 
  outr[5] += incr[5]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += -1.0*incr[3]*dv11; 
  outl[4] += incr[4]*dv11; 
  outl[5] += -1.0*incr[5]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[6]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 

  double abar2[6]; 


  abar1[0] = E1[0]+wv3*B0[0]-wv1*B2[0]; 
  abar1[1] = E1[1]+wv3*B0[1]-wv1*B2[1]; 
  abar1[2] = E1[2]+wv3*B0[2]-wv1*B2[2]; 
  abar1[3] = E1[3]+wv3*B0[3]-wv1*B2[3]; 
  abar1[4] = E1[4]+wv3*B0[4]-wv1*B2[4]; 
  abar1[5] = E1[5]+wv3*B0[5]-wv1*B2[5]; 

  double incr[21]; 

  double favg[21]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = -1*fr[4]+fl[4]; 
  favg[5] = 1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = -1*fr[10]+fl[10]; 
  favg[11] = -1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = -1*fr[4]-fl[4]; 
  fjump[5] = 1*fr[5]-fl[5]; 
  fjump[6] = 1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = -1*fr[9]-fl[9]; 
  fjump[10] = -1*fr[10]-fl[10]; 
  fjump[11] = -1*fr[11]-fl[11]; 
  fjump[12] = 1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = 1*fr[14]-fl[14]; 
  fjump[15] = -1*fr[15]-fl[15]; 
  fjump[16] = 1*fr[16]-fl[16]; 
  fjump[17] = 1*fr[17]-fl[17]; 
  fjump[18] = 1*fr[18]-fl[18]; 
  fjump[19] = 1*fr[19]-fl[19]; 
  fjump[20] = 1*fr[20]-fl[20]; 
  const double amid = (-0.5590169943749475*abar1[5])-0.5590169943749475*abar1[4]+0.5*abar1[0]; 
  incr[0] = 0.0625*B0[0]*favg[15]*dv3+0.03608439182435161*B0[2]*favg[13]*dv3+0.03608439182435161*B0[1]*favg[12]*dv3+0.03608439182435161*B0[0]*favg[5]*dv3-0.0625*B2[0]*favg[11]*dv1-0.03608439182435161*B2[2]*favg[8]*dv1-0.03608439182435161*B2[1]*favg[7]*dv1-0.03608439182435161*B2[0]*favg[3]*dv1-0.5590169943749475*fjump[19]*amax-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.2795084971874737*abar1[0]*favg[19]+0.125*abar1[5]*favg[17]+0.125*abar1[4]*favg[16]+0.2165063509461096*abar1[2]*favg[10]+0.2165063509461096*abar1[1]*favg[9]+0.125*abar1[3]*favg[6]+0.2165063509461096*abar1[0]*favg[4]+0.125*abar1[2]*favg[2]+0.125*abar1[1]*favg[1]+0.125*abar1[0]*favg[0]; 
  incr[1] = 0.0625*B0[1]*favg[15]*dv3+0.03608439182435161*B0[3]*favg[13]*dv3+0.03227486121839514*B0[4]*favg[12]*dv3+0.03608439182435161*B0[0]*favg[12]*dv3+0.03608439182435161*B0[1]*favg[5]*dv3-0.0625*B2[1]*favg[11]*dv1-0.03608439182435161*B2[3]*favg[8]*dv1-0.03227486121839514*B2[4]*favg[7]*dv1-0.03608439182435161*B2[0]*favg[7]*dv1-0.03608439182435161*B2[1]*favg[3]*dv1-0.4330127018922193*fjump[9]*amax-0.25*fjump[1]*amax+0.2795084971874738*abar1[1]*favg[19]+0.1118033988749895*abar1[1]*favg[16]+0.2165063509461096*abar1[3]*favg[10]+0.1936491673103708*abar1[4]*favg[9]+0.2165063509461096*abar1[0]*favg[9]+0.125*abar1[2]*favg[6]+0.2165063509461096*abar1[1]*favg[4]+0.1118033988749895*favg[1]*abar1[4]+0.125*favg[2]*abar1[3]+0.125*abar1[0]*favg[1]+0.125*favg[0]*abar1[1]; 
  incr[2] = 0.0625*B0[2]*favg[15]*dv3+0.03227486121839514*B0[5]*favg[13]*dv3+0.03608439182435161*B0[0]*favg[13]*dv3+0.03608439182435161*B0[3]*favg[12]*dv3+0.03608439182435161*B0[2]*favg[5]*dv3-0.0625*B2[2]*favg[11]*dv1-0.03227486121839514*B2[5]*favg[8]*dv1-0.03608439182435161*B2[0]*favg[8]*dv1-0.03608439182435161*B2[3]*favg[7]*dv1-0.03608439182435161*B2[2]*favg[3]*dv1-0.4330127018922193*fjump[10]*amax-0.25*fjump[2]*amax+0.2795084971874738*abar1[2]*favg[19]+0.1118033988749895*abar1[2]*favg[17]+0.1936491673103708*abar1[5]*favg[10]+0.2165063509461096*abar1[0]*favg[10]+0.2165063509461096*abar1[3]*favg[9]+0.125*abar1[1]*favg[6]+0.1118033988749895*favg[2]*abar1[5]+0.2165063509461096*abar1[2]*favg[4]+0.125*favg[1]*abar1[3]+0.125*abar1[0]*favg[2]+0.125*favg[0]*abar1[2]; 
  incr[3] = 0.03608439182435161*B0[0]*favg[14]*dv3-0.08068715304598786*B2[0]*favg[19]*dv1-0.03227486121839514*B2[0]*favg[18]*dv1-0.03608439182435161*B2[5]*favg[17]*dv1-0.03608439182435161*B2[4]*favg[16]*dv1-0.0625*B2[2]*favg[10]*dv1-0.0625*B2[1]*favg[9]*dv1-0.03608439182435161*B2[3]*favg[6]*dv1-0.0625*B2[0]*favg[4]*dv1-0.03608439182435161*favg[2]*B2[2]*dv1-0.03608439182435161*favg[1]*B2[1]*dv1-0.03608439182435161*favg[0]*B2[0]*dv1-0.4330127018922193*fjump[11]*amax-0.25*fjump[3]*amax+0.2165063509461096*abar1[0]*favg[11]+0.125*abar1[2]*favg[8]+0.125*abar1[1]*favg[7]+0.125*abar1[0]*favg[3]; 
  incr[4] = (-0.1082531754730548*B0[0]*favg[15]*dv3)-0.0625*B0[2]*favg[13]*dv3-0.0625*B0[1]*favg[12]*dv3-0.0625*B0[0]*favg[5]*dv3+0.1082531754730548*B2[0]*favg[11]*dv1+0.0625*B2[2]*favg[8]*dv1+0.0625*B2[1]*favg[7]*dv1+0.0625*B2[0]*favg[3]*dv1+0.9682458365518543*fjump[19]*amax+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.4841229182759271*abar1[0]*favg[19]-0.2165063509461096*abar1[5]*favg[17]-0.2165063509461096*abar1[4]*favg[16]-0.375*abar1[2]*favg[10]-0.375*abar1[1]*favg[9]-0.2165063509461096*abar1[3]*favg[6]-0.375*abar1[0]*favg[4]-0.2165063509461096*abar1[2]*favg[2]-0.2165063509461096*abar1[1]*favg[1]-0.2165063509461096*abar1[0]*favg[0]; 
  incr[5] = 0.03227486121839514*B0[0]*favg[20]*dv3+0.08068715304598786*B0[0]*favg[19]*dv3+0.03608439182435161*B0[5]*favg[17]*dv3+0.03608439182435161*B0[4]*favg[16]*dv3+0.0625*B0[2]*favg[10]*dv3+0.0625*B0[1]*favg[9]*dv3+0.03608439182435161*B0[3]*favg[6]*dv3+0.0625*B0[0]*favg[4]*dv3+0.03608439182435161*favg[2]*B0[2]*dv3+0.03608439182435161*favg[1]*B0[1]*dv3+0.03608439182435161*favg[0]*B0[0]*dv3-0.03608439182435161*B2[0]*favg[14]*dv1-0.4330127018922193*fjump[15]*amax-0.25*fjump[5]*amax+0.2165063509461096*abar1[0]*favg[15]+0.125*abar1[2]*favg[13]+0.125*abar1[1]*favg[12]+0.125*abar1[0]*favg[5]; 
  incr[6] = 0.0625*B0[3]*favg[15]*dv3+0.03608439182435161*B0[1]*favg[13]*dv3+0.03608439182435161*B0[2]*favg[12]*dv3+0.03608439182435161*B0[3]*favg[5]*dv3-0.0625*B2[3]*favg[11]*dv1-0.03608439182435161*B2[1]*favg[8]*dv1-0.03608439182435161*B2[2]*favg[7]*dv1-0.03608439182435161*favg[3]*B2[3]*dv1-0.25*fjump[6]*amax+0.2795084971874738*abar1[3]*favg[19]+0.1118033988749895*abar1[3]*favg[17]+0.1118033988749895*abar1[3]*favg[16]+0.2165063509461096*abar1[1]*favg[10]+0.2165063509461096*abar1[2]*favg[9]+0.1118033988749895*abar1[5]*favg[6]+0.1118033988749895*abar1[4]*favg[6]+0.125*abar1[0]*favg[6]+0.2165063509461096*abar1[3]*favg[4]+0.125*favg[0]*abar1[3]+0.125*abar1[1]*favg[2]+0.125*favg[1]*abar1[2]; 
  incr[7] = 0.03608439182435161*B0[1]*favg[14]*dv3-0.08068715304598786*B2[1]*favg[19]*dv1-0.03227486121839514*B2[1]*favg[18]*dv1-0.03227486121839514*B2[1]*favg[16]*dv1-0.0625*B2[3]*favg[10]*dv1-0.05590169943749475*B2[4]*favg[9]*dv1-0.0625*B2[0]*favg[9]*dv1-0.03608439182435161*B2[2]*favg[6]*dv1-0.03227486121839514*favg[1]*B2[4]*dv1-0.0625*B2[1]*favg[4]*dv1-0.03608439182435161*favg[2]*B2[3]*dv1-0.03608439182435161*favg[0]*B2[1]*dv1-0.03608439182435161*B2[0]*favg[1]*dv1-0.25*fjump[7]*amax+0.2165063509461096*abar1[1]*favg[11]+0.125*abar1[3]*favg[8]+0.1118033988749895*abar1[4]*favg[7]+0.125*abar1[0]*favg[7]+0.125*abar1[1]*favg[3]; 
  incr[8] = 0.03608439182435161*B0[2]*favg[14]*dv3-0.08068715304598786*B2[2]*favg[19]*dv1-0.03227486121839514*B2[2]*favg[18]*dv1-0.03227486121839514*B2[2]*favg[17]*dv1-0.05590169943749475*B2[5]*favg[10]*dv1-0.0625*B2[0]*favg[10]*dv1-0.0625*B2[3]*favg[9]*dv1-0.03608439182435161*B2[1]*favg[6]*dv1-0.03227486121839514*favg[2]*B2[5]*dv1-0.0625*B2[2]*favg[4]*dv1-0.03608439182435161*favg[1]*B2[3]*dv1-0.03608439182435161*favg[0]*B2[2]*dv1-0.03608439182435161*B2[0]*favg[2]*dv1-0.25*fjump[8]*amax+0.2165063509461096*abar1[2]*favg[11]+0.1118033988749895*abar1[5]*favg[8]+0.125*abar1[0]*favg[8]+0.125*abar1[3]*favg[7]+0.125*abar1[2]*favg[3]; 
  incr[9] = (-0.1082531754730548*B0[1]*favg[15]*dv3)-0.0625*B0[3]*favg[13]*dv3-0.05590169943749475*B0[4]*favg[12]*dv3-0.0625*B0[0]*favg[12]*dv3-0.0625*B0[1]*favg[5]*dv3+0.1082531754730548*B2[1]*favg[11]*dv1+0.0625*B2[3]*favg[8]*dv1+0.05590169943749475*B2[4]*favg[7]*dv1+0.0625*B2[0]*favg[7]*dv1+0.0625*B2[1]*favg[3]*dv1+0.75*fjump[9]*amax+0.4330127018922193*fjump[1]*amax-0.4841229182759272*abar1[1]*favg[19]-0.1936491673103709*abar1[1]*favg[16]-0.375*abar1[3]*favg[10]-0.3354101966249685*abar1[4]*favg[9]-0.375*abar1[0]*favg[9]-0.2165063509461096*abar1[2]*favg[6]-0.375*abar1[1]*favg[4]-0.1936491673103709*favg[1]*abar1[4]-0.2165063509461096*favg[2]*abar1[3]-0.2165063509461096*abar1[0]*favg[1]-0.2165063509461096*favg[0]*abar1[1]; 
  incr[10] = (-0.1082531754730548*B0[2]*favg[15]*dv3)-0.05590169943749475*B0[5]*favg[13]*dv3-0.0625*B0[0]*favg[13]*dv3-0.0625*B0[3]*favg[12]*dv3-0.0625*B0[2]*favg[5]*dv3+0.1082531754730548*B2[2]*favg[11]*dv1+0.05590169943749475*B2[5]*favg[8]*dv1+0.0625*B2[0]*favg[8]*dv1+0.0625*B2[3]*favg[7]*dv1+0.0625*B2[2]*favg[3]*dv1+0.75*fjump[10]*amax+0.4330127018922193*fjump[2]*amax-0.4841229182759272*abar1[2]*favg[19]-0.1936491673103709*abar1[2]*favg[17]-0.3354101966249685*abar1[5]*favg[10]-0.375*abar1[0]*favg[10]-0.375*abar1[3]*favg[9]-0.2165063509461096*abar1[1]*favg[6]-0.1936491673103709*favg[2]*abar1[5]-0.375*abar1[2]*favg[4]-0.2165063509461096*favg[1]*abar1[3]-0.2165063509461096*abar1[0]*favg[2]-0.2165063509461096*favg[0]*abar1[2]; 
  incr[11] = (-0.0625*B0[0]*favg[14]*dv3)+0.1397542485937369*B2[0]*favg[19]*dv1+0.05590169943749475*B2[0]*favg[18]*dv1+0.0625*B2[5]*favg[17]*dv1+0.0625*B2[4]*favg[16]*dv1+0.1082531754730548*B2[2]*favg[10]*dv1+0.1082531754730548*B2[1]*favg[9]*dv1+0.0625*B2[3]*favg[6]*dv1+0.1082531754730548*B2[0]*favg[4]*dv1+0.0625*favg[2]*B2[2]*dv1+0.0625*favg[1]*B2[1]*dv1+0.0625*favg[0]*B2[0]*dv1+0.75*fjump[11]*amax+0.4330127018922193*fjump[3]*amax-0.375*abar1[0]*favg[11]-0.2165063509461096*abar1[2]*favg[8]-0.2165063509461096*abar1[1]*favg[7]-0.2165063509461096*abar1[0]*favg[3]; 
  incr[12] = 0.03227486121839514*B0[1]*favg[20]*dv3+0.08068715304598786*B0[1]*favg[19]*dv3+0.03227486121839514*B0[1]*favg[16]*dv3+0.0625*B0[3]*favg[10]*dv3+0.05590169943749475*B0[4]*favg[9]*dv3+0.0625*B0[0]*favg[9]*dv3+0.03608439182435161*B0[2]*favg[6]*dv3+0.03227486121839514*favg[1]*B0[4]*dv3+0.0625*B0[1]*favg[4]*dv3+0.03608439182435161*favg[2]*B0[3]*dv3+0.03608439182435161*favg[0]*B0[1]*dv3+0.03608439182435161*B0[0]*favg[1]*dv3-0.03608439182435161*B2[1]*favg[14]*dv1-0.25*fjump[12]*amax+0.2165063509461096*abar1[1]*favg[15]+0.125*abar1[3]*favg[13]+0.1118033988749895*abar1[4]*favg[12]+0.125*abar1[0]*favg[12]+0.125*abar1[1]*favg[5]; 
  incr[13] = 0.03227486121839514*B0[2]*favg[20]*dv3+0.08068715304598786*B0[2]*favg[19]*dv3+0.03227486121839514*B0[2]*favg[17]*dv3+0.05590169943749475*B0[5]*favg[10]*dv3+0.0625*B0[0]*favg[10]*dv3+0.0625*B0[3]*favg[9]*dv3+0.03608439182435161*B0[1]*favg[6]*dv3+0.03227486121839514*favg[2]*B0[5]*dv3+0.0625*B0[2]*favg[4]*dv3+0.03608439182435161*favg[1]*B0[3]*dv3+0.03608439182435161*favg[0]*B0[2]*dv3+0.03608439182435161*B0[0]*favg[2]*dv3-0.03608439182435161*B2[2]*favg[14]*dv1-0.25*fjump[13]*amax+0.2165063509461096*abar1[2]*favg[15]+0.1118033988749895*abar1[5]*favg[13]+0.125*abar1[0]*favg[13]+0.125*abar1[3]*favg[12]+0.125*abar1[2]*favg[5]; 
  incr[14] = 0.0625*B0[0]*favg[11]*dv3+0.03608439182435161*B0[2]*favg[8]*dv3+0.03608439182435161*B0[1]*favg[7]*dv3+0.03608439182435161*B0[0]*favg[3]*dv3-0.0625*B2[0]*favg[15]*dv1-0.03608439182435161*B2[2]*favg[13]*dv1-0.03608439182435161*B2[1]*favg[12]*dv1-0.03608439182435161*B2[0]*favg[5]*dv1-0.25*fjump[14]*amax+0.125*abar1[0]*favg[14]; 
  incr[15] = (-0.05590169943749475*B0[0]*favg[20]*dv3)-0.1397542485937369*B0[0]*favg[19]*dv3-0.0625*B0[5]*favg[17]*dv3-0.0625*B0[4]*favg[16]*dv3-0.1082531754730548*B0[2]*favg[10]*dv3-0.1082531754730548*B0[1]*favg[9]*dv3-0.0625*B0[3]*favg[6]*dv3-0.1082531754730548*B0[0]*favg[4]*dv3-0.0625*favg[2]*B0[2]*dv3-0.0625*favg[1]*B0[1]*dv3-0.0625*favg[0]*B0[0]*dv3+0.0625*B2[0]*favg[14]*dv1+0.75*fjump[15]*amax+0.4330127018922193*fjump[5]*amax-0.375*abar1[0]*favg[15]-0.2165063509461096*abar1[2]*favg[13]-0.2165063509461096*abar1[1]*favg[12]-0.2165063509461096*abar1[0]*favg[5]; 
  incr[16] = 0.0625*B0[4]*favg[15]*dv3+0.03227486121839514*B0[1]*favg[12]*dv3+0.03608439182435162*B0[4]*favg[5]*dv3-0.0625*B2[4]*favg[11]*dv1-0.03227486121839514*B2[1]*favg[7]*dv1-0.03608439182435162*favg[3]*B2[4]*dv1-0.25*fjump[16]*amax+0.2795084971874738*abar1[4]*favg[19]+0.0798595706249925*abar1[4]*favg[16]+0.125*abar1[0]*favg[16]+0.1936491673103709*abar1[1]*favg[9]+0.1118033988749895*abar1[3]*favg[6]+0.2165063509461096*abar1[4]*favg[4]+0.125*favg[0]*abar1[4]+0.1118033988749895*abar1[1]*favg[1]; 
  incr[17] = 0.0625*B0[5]*favg[15]*dv3+0.03227486121839514*B0[2]*favg[13]*dv3+0.03608439182435162*favg[5]*B0[5]*dv3-0.0625*B2[5]*favg[11]*dv1-0.03227486121839514*B2[2]*favg[8]*dv1-0.03608439182435162*favg[3]*B2[5]*dv1-0.25*fjump[17]*amax+0.2795084971874738*abar1[5]*favg[19]+0.0798595706249925*abar1[5]*favg[17]+0.125*abar1[0]*favg[17]+0.1936491673103709*abar1[2]*favg[10]+0.1118033988749895*abar1[3]*favg[6]+0.2165063509461096*favg[4]*abar1[5]+0.125*favg[0]*abar1[5]+0.1118033988749895*abar1[2]*favg[2]; 
  incr[18] = (-0.05590169943749475*B2[0]*favg[11]*dv1)-0.03227486121839514*B2[2]*favg[8]*dv1-0.03227486121839514*B2[1]*favg[7]*dv1-0.03227486121839514*B2[0]*favg[3]*dv1-0.25*fjump[18]*amax+0.125*abar1[0]*favg[18]; 
  incr[19] = 0.1397542485937369*B0[0]*favg[15]*dv3+0.08068715304598785*B0[2]*favg[13]*dv3+0.08068715304598785*B0[1]*favg[12]*dv3+0.08068715304598785*B0[0]*favg[5]*dv3-0.1397542485937369*B2[0]*favg[11]*dv1-0.08068715304598785*B2[2]*favg[8]*dv1-0.08068715304598785*B2[1]*favg[7]*dv1-0.08068715304598785*B2[0]*favg[3]*dv1-1.25*fjump[19]*amax-0.9682458365518543*fjump[4]*amax-0.5590169943749475*fjump[0]*amax+0.625*abar1[0]*favg[19]+0.2795084971874737*abar1[5]*favg[17]+0.2795084971874737*abar1[4]*favg[16]+0.484122918275927*abar1[2]*favg[10]+0.484122918275927*abar1[1]*favg[9]+0.2795084971874737*abar1[3]*favg[6]+0.484122918275927*abar1[0]*favg[4]+0.2795084971874737*abar1[2]*favg[2]+0.2795084971874737*abar1[1]*favg[1]+0.2795084971874737*abar1[0]*favg[0]; 
  incr[20] = 0.05590169943749475*B0[0]*favg[15]*dv3+0.03227486121839514*B0[2]*favg[13]*dv3+0.03227486121839514*B0[1]*favg[12]*dv3+0.03227486121839514*B0[0]*favg[5]*dv3-0.25*fjump[20]*amax+0.125*abar1[0]*favg[20]; 

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
  outr[15] += incr[15]*dv11; 
  outr[16] += incr[16]*dv11; 
  outr[17] += incr[17]*dv11; 
  outr[18] += incr[18]*dv11; 
  outr[19] += incr[19]*dv11; 
  outr[20] += incr[20]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += -1.0*incr[3]*dv11; 
  outl[4] += incr[4]*dv11; 
  outl[5] += -1.0*incr[5]*dv11; 
  outl[6] += -1.0*incr[6]*dv11; 
  outl[7] += -1.0*incr[7]*dv11; 
  outl[8] += -1.0*incr[8]*dv11; 
  outl[9] += incr[9]*dv11; 
  outl[10] += incr[10]*dv11; 
  outl[11] += incr[11]*dv11; 
  outl[12] += -1.0*incr[12]*dv11; 
  outl[13] += -1.0*incr[13]*dv11; 
  outl[14] += -1.0*incr[14]*dv11; 
  outl[15] += incr[15]*dv11; 
  outl[16] += -1.0*incr[16]*dv11; 
  outl[17] += -1.0*incr[17]*dv11; 
  outl[18] += -1.0*incr[18]*dv11; 
  outl[19] += -1.0*incr[19]*dv11; 
  outl[20] += -1.0*incr[20]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VZ_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12 = 2/dxv[4]; 
  const double *E2 = &EM[6]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar2[0] = E2[0]+wv1*B1[0]-wv2*B0[0]; 
  abar2[1] = E2[1]+wv1*B1[1]-wv2*B0[1]; 
  abar2[2] = E2[2]+wv1*B1[2]-wv2*B0[2]; 

  double incr[6]; 

  double favg[6]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  double fjump[6]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = -1*fr[5]-fl[5]; 
  const double amid = 0.5*abar2[0]; 
  incr[0] = (-0.03608439182435161*B0[0]*favg[4]*dv2)+0.03608439182435161*B1[0]*favg[3]*dv1-0.4330127018922193*fjump[5]*amax-0.25*fjump[0]*amax+0.2165063509461096*abar2[0]*favg[5]+0.125*abar2[2]*favg[2]+0.125*abar2[1]*favg[1]+0.125*abar2[0]*favg[0]; 
  incr[1] = (-0.03608439182435161*B0[1]*favg[4]*dv2)+0.03608439182435161*B1[1]*favg[3]*dv1-0.25*fjump[1]*amax+0.2165063509461096*abar2[1]*favg[5]+0.125*abar2[0]*favg[1]+0.125*favg[0]*abar2[1]; 
  incr[2] = (-0.03608439182435161*B0[2]*favg[4]*dv2)+0.03608439182435161*B1[2]*favg[3]*dv1-0.25*fjump[2]*amax+0.2165063509461096*abar2[2]*favg[5]+0.125*abar2[0]*favg[2]+0.125*favg[0]*abar2[2]; 
  incr[3] = 0.0625*B1[0]*favg[5]*dv1+0.03608439182435161*favg[2]*B1[2]*dv1+0.03608439182435161*favg[1]*B1[1]*dv1+0.03608439182435161*favg[0]*B1[0]*dv1-0.25*fjump[3]*amax+0.125*abar2[0]*favg[3]; 
  incr[4] = (-0.0625*B0[0]*favg[5]*dv2)-0.03608439182435161*favg[2]*B0[2]*dv2-0.03608439182435161*favg[1]*B0[1]*dv2-0.03608439182435161*favg[0]*B0[0]*dv2-0.25*fjump[4]*amax+0.125*abar2[0]*favg[4]; 
  incr[5] = 0.0625*B0[0]*favg[4]*dv2-0.0625*B1[0]*favg[3]*dv1+0.75*fjump[5]*amax+0.4330127018922193*fjump[0]*amax-0.375*abar2[0]*favg[5]-0.2165063509461096*abar2[2]*favg[2]-0.2165063509461096*abar2[1]*favg[1]-0.2165063509461096*abar2[0]*favg[0]; 

  outr[0] += incr[0]*dv12; 
  outr[1] += incr[1]*dv12; 
  outr[2] += incr[2]*dv12; 
  outr[3] += incr[3]*dv12; 
  outr[4] += incr[4]*dv12; 
  outr[5] += incr[5]*dv12; 

  outl[0] += -1.0*incr[0]*dv12; 
  outl[1] += -1.0*incr[1]*dv12; 
  outl[2] += -1.0*incr[2]*dv12; 
  outl[3] += -1.0*incr[3]*dv12; 
  outl[4] += -1.0*incr[4]*dv12; 
  outl[5] += incr[5]*dv12; 
return std::abs(amid); 
} 
double VlasovSurfElcMag2x3vMax_VZ_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12 = 2/dxv[4]; 
  const double *E2 = &EM[12]; 

  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 

  double abar2[6]; 


  abar2[0] = E2[0]+wv1*B1[0]-wv2*B0[0]; 
  abar2[1] = E2[1]+wv1*B1[1]-wv2*B0[1]; 
  abar2[2] = E2[2]+wv1*B1[2]-wv2*B0[2]; 
  abar2[3] = E2[3]+wv1*B1[3]-wv2*B0[3]; 
  abar2[4] = E2[4]+wv1*B1[4]-wv2*B0[4]; 
  abar2[5] = E2[5]+wv1*B1[5]-wv2*B0[5]; 

  double incr[21]; 

  double favg[21]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = 1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = 1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = 1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = -1*fr[12]+fl[12]; 
  favg[13] = -1*fr[13]+fl[13]; 
  favg[14] = -1*fr[14]+fl[14]; 
  favg[15] = -1*fr[15]+fl[15]; 
  favg[16] = 1*fr[16]+fl[16]; 
  favg[17] = 1*fr[17]+fl[17]; 
  favg[18] = 1*fr[18]+fl[18]; 
  favg[19] = 1*fr[19]+fl[19]; 
  favg[20] = 1*fr[20]+fl[20]; 
  double fjump[21]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = 1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = -1*fr[5]-fl[5]; 
  fjump[6] = 1*fr[6]-fl[6]; 
  fjump[7] = 1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = 1*fr[9]-fl[9]; 
  fjump[10] = 1*fr[10]-fl[10]; 
  fjump[11] = 1*fr[11]-fl[11]; 
  fjump[12] = -1*fr[12]-fl[12]; 
  fjump[13] = -1*fr[13]-fl[13]; 
  fjump[14] = -1*fr[14]-fl[14]; 
  fjump[15] = -1*fr[15]-fl[15]; 
  fjump[16] = 1*fr[16]-fl[16]; 
  fjump[17] = 1*fr[17]-fl[17]; 
  fjump[18] = 1*fr[18]-fl[18]; 
  fjump[19] = 1*fr[19]-fl[19]; 
  fjump[20] = 1*fr[20]-fl[20]; 
  const double amid = (-0.5590169943749475*abar2[5])-0.5590169943749475*abar2[4]+0.5*abar2[0]; 
  incr[0] = (-0.0625*B0[0]*favg[15]*dv2)-0.03608439182435161*B0[2]*favg[10]*dv2-0.03608439182435161*B0[1]*favg[9]*dv2-0.03608439182435161*B0[0]*favg[4]*dv2+0.0625*B1[0]*favg[14]*dv1+0.03608439182435161*B1[2]*favg[8]*dv1+0.03608439182435161*B1[1]*favg[7]*dv1+0.03608439182435161*B1[0]*favg[3]*dv1-0.5590169943749475*fjump[20]*amax-0.4330127018922193*fjump[5]*amax-0.25*fjump[0]*amax+0.2795084971874737*abar2[0]*favg[20]+0.125*abar2[5]*favg[17]+0.125*abar2[4]*favg[16]+0.2165063509461096*abar2[2]*favg[13]+0.2165063509461096*abar2[1]*favg[12]+0.125*abar2[3]*favg[6]+0.2165063509461096*abar2[0]*favg[5]+0.125*abar2[2]*favg[2]+0.125*abar2[1]*favg[1]+0.125*abar2[0]*favg[0]; 
  incr[1] = (-0.0625*B0[1]*favg[15]*dv2)-0.03608439182435161*B0[3]*favg[10]*dv2-0.03227486121839514*B0[4]*favg[9]*dv2-0.03608439182435161*B0[0]*favg[9]*dv2-0.03608439182435161*B0[1]*favg[4]*dv2+0.0625*B1[1]*favg[14]*dv1+0.03608439182435161*B1[3]*favg[8]*dv1+0.03227486121839514*B1[4]*favg[7]*dv1+0.03608439182435161*B1[0]*favg[7]*dv1+0.03608439182435161*B1[1]*favg[3]*dv1-0.4330127018922193*fjump[12]*amax-0.25*fjump[1]*amax+0.2795084971874738*abar2[1]*favg[20]+0.1118033988749895*abar2[1]*favg[16]+0.2165063509461096*abar2[3]*favg[13]+0.1936491673103708*abar2[4]*favg[12]+0.2165063509461096*abar2[0]*favg[12]+0.125*abar2[2]*favg[6]+0.2165063509461096*abar2[1]*favg[5]+0.1118033988749895*favg[1]*abar2[4]+0.125*favg[2]*abar2[3]+0.125*abar2[0]*favg[1]+0.125*favg[0]*abar2[1]; 
  incr[2] = (-0.0625*B0[2]*favg[15]*dv2)-0.03227486121839514*B0[5]*favg[10]*dv2-0.03608439182435161*B0[0]*favg[10]*dv2-0.03608439182435161*B0[3]*favg[9]*dv2-0.03608439182435161*B0[2]*favg[4]*dv2+0.0625*B1[2]*favg[14]*dv1+0.03227486121839514*B1[5]*favg[8]*dv1+0.03608439182435161*B1[0]*favg[8]*dv1+0.03608439182435161*B1[3]*favg[7]*dv1+0.03608439182435161*B1[2]*favg[3]*dv1-0.4330127018922193*fjump[13]*amax-0.25*fjump[2]*amax+0.2795084971874738*abar2[2]*favg[20]+0.1118033988749895*abar2[2]*favg[17]+0.1936491673103708*abar2[5]*favg[13]+0.2165063509461096*abar2[0]*favg[13]+0.2165063509461096*abar2[3]*favg[12]+0.125*abar2[1]*favg[6]+0.2165063509461096*abar2[2]*favg[5]+0.1118033988749895*favg[2]*abar2[5]+0.125*favg[1]*abar2[3]+0.125*abar2[0]*favg[2]+0.125*favg[0]*abar2[2]; 
  incr[3] = (-0.03608439182435161*B0[0]*favg[11]*dv2)+0.08068715304598786*B1[0]*favg[20]*dv1+0.03227486121839514*B1[0]*favg[18]*dv1+0.03608439182435161*B1[5]*favg[17]*dv1+0.03608439182435161*B1[4]*favg[16]*dv1+0.0625*B1[2]*favg[13]*dv1+0.0625*B1[1]*favg[12]*dv1+0.03608439182435161*B1[3]*favg[6]*dv1+0.0625*B1[0]*favg[5]*dv1+0.03608439182435161*favg[2]*B1[2]*dv1+0.03608439182435161*favg[1]*B1[1]*dv1+0.03608439182435161*favg[0]*B1[0]*dv1-0.4330127018922193*fjump[14]*amax-0.25*fjump[3]*amax+0.2165063509461096*abar2[0]*favg[14]+0.125*abar2[2]*favg[8]+0.125*abar2[1]*favg[7]+0.125*abar2[0]*favg[3]; 
  incr[4] = (-0.08068715304598786*B0[0]*favg[20]*dv2)-0.03227486121839514*B0[0]*favg[19]*dv2-0.03608439182435161*B0[5]*favg[17]*dv2-0.03608439182435161*B0[4]*favg[16]*dv2-0.0625*B0[2]*favg[13]*dv2-0.0625*B0[1]*favg[12]*dv2-0.03608439182435161*B0[3]*favg[6]*dv2-0.0625*B0[0]*favg[5]*dv2-0.03608439182435161*favg[2]*B0[2]*dv2-0.03608439182435161*favg[1]*B0[1]*dv2-0.03608439182435161*favg[0]*B0[0]*dv2+0.03608439182435161*B1[0]*favg[11]*dv1-0.4330127018922193*fjump[15]*amax-0.25*fjump[4]*amax+0.2165063509461096*abar2[0]*favg[15]+0.125*abar2[2]*favg[10]+0.125*abar2[1]*favg[9]+0.125*abar2[0]*favg[4]; 
  incr[5] = 0.1082531754730548*B0[0]*favg[15]*dv2+0.0625*B0[2]*favg[10]*dv2+0.0625*B0[1]*favg[9]*dv2+0.0625*B0[0]*favg[4]*dv2-0.1082531754730548*B1[0]*favg[14]*dv1-0.0625*B1[2]*favg[8]*dv1-0.0625*B1[1]*favg[7]*dv1-0.0625*B1[0]*favg[3]*dv1+0.9682458365518543*fjump[20]*amax+0.75*fjump[5]*amax+0.4330127018922193*fjump[0]*amax-0.4841229182759271*abar2[0]*favg[20]-0.2165063509461096*abar2[5]*favg[17]-0.2165063509461096*abar2[4]*favg[16]-0.375*abar2[2]*favg[13]-0.375*abar2[1]*favg[12]-0.2165063509461096*abar2[3]*favg[6]-0.375*abar2[0]*favg[5]-0.2165063509461096*abar2[2]*favg[2]-0.2165063509461096*abar2[1]*favg[1]-0.2165063509461096*abar2[0]*favg[0]; 
  incr[6] = (-0.0625*B0[3]*favg[15]*dv2)-0.03608439182435161*B0[1]*favg[10]*dv2-0.03608439182435161*B0[2]*favg[9]*dv2-0.03608439182435161*B0[3]*favg[4]*dv2+0.0625*B1[3]*favg[14]*dv1+0.03608439182435161*B1[1]*favg[8]*dv1+0.03608439182435161*B1[2]*favg[7]*dv1+0.03608439182435161*favg[3]*B1[3]*dv1-0.25*fjump[6]*amax+0.2795084971874738*abar2[3]*favg[20]+0.1118033988749895*abar2[3]*favg[17]+0.1118033988749895*abar2[3]*favg[16]+0.2165063509461096*abar2[1]*favg[13]+0.2165063509461096*abar2[2]*favg[12]+0.1118033988749895*abar2[5]*favg[6]+0.1118033988749895*abar2[4]*favg[6]+0.125*abar2[0]*favg[6]+0.2165063509461096*abar2[3]*favg[5]+0.125*favg[0]*abar2[3]+0.125*abar2[1]*favg[2]+0.125*favg[1]*abar2[2]; 
  incr[7] = (-0.03608439182435161*B0[1]*favg[11]*dv2)+0.08068715304598786*B1[1]*favg[20]*dv1+0.03227486121839514*B1[1]*favg[18]*dv1+0.03227486121839514*B1[1]*favg[16]*dv1+0.0625*B1[3]*favg[13]*dv1+0.05590169943749475*B1[4]*favg[12]*dv1+0.0625*B1[0]*favg[12]*dv1+0.03608439182435161*B1[2]*favg[6]*dv1+0.0625*B1[1]*favg[5]*dv1+0.03227486121839514*favg[1]*B1[4]*dv1+0.03608439182435161*favg[2]*B1[3]*dv1+0.03608439182435161*favg[0]*B1[1]*dv1+0.03608439182435161*B1[0]*favg[1]*dv1-0.25*fjump[7]*amax+0.2165063509461096*abar2[1]*favg[14]+0.125*abar2[3]*favg[8]+0.1118033988749895*abar2[4]*favg[7]+0.125*abar2[0]*favg[7]+0.125*abar2[1]*favg[3]; 
  incr[8] = (-0.03608439182435161*B0[2]*favg[11]*dv2)+0.08068715304598786*B1[2]*favg[20]*dv1+0.03227486121839514*B1[2]*favg[18]*dv1+0.03227486121839514*B1[2]*favg[17]*dv1+0.05590169943749475*B1[5]*favg[13]*dv1+0.0625*B1[0]*favg[13]*dv1+0.0625*B1[3]*favg[12]*dv1+0.03608439182435161*B1[1]*favg[6]*dv1+0.03227486121839514*favg[2]*B1[5]*dv1+0.0625*B1[2]*favg[5]*dv1+0.03608439182435161*favg[1]*B1[3]*dv1+0.03608439182435161*favg[0]*B1[2]*dv1+0.03608439182435161*B1[0]*favg[2]*dv1-0.25*fjump[8]*amax+0.2165063509461096*abar2[2]*favg[14]+0.1118033988749895*abar2[5]*favg[8]+0.125*abar2[0]*favg[8]+0.125*abar2[3]*favg[7]+0.125*abar2[2]*favg[3]; 
  incr[9] = (-0.08068715304598786*B0[1]*favg[20]*dv2)-0.03227486121839514*B0[1]*favg[19]*dv2-0.03227486121839514*B0[1]*favg[16]*dv2-0.0625*B0[3]*favg[13]*dv2-0.05590169943749475*B0[4]*favg[12]*dv2-0.0625*B0[0]*favg[12]*dv2-0.03608439182435161*B0[2]*favg[6]*dv2-0.0625*B0[1]*favg[5]*dv2-0.03227486121839514*favg[1]*B0[4]*dv2-0.03608439182435161*favg[2]*B0[3]*dv2-0.03608439182435161*favg[0]*B0[1]*dv2-0.03608439182435161*B0[0]*favg[1]*dv2+0.03608439182435161*B1[1]*favg[11]*dv1-0.25*fjump[9]*amax+0.2165063509461096*abar2[1]*favg[15]+0.125*abar2[3]*favg[10]+0.1118033988749895*abar2[4]*favg[9]+0.125*abar2[0]*favg[9]+0.125*abar2[1]*favg[4]; 
  incr[10] = (-0.08068715304598786*B0[2]*favg[20]*dv2)-0.03227486121839514*B0[2]*favg[19]*dv2-0.03227486121839514*B0[2]*favg[17]*dv2-0.05590169943749475*B0[5]*favg[13]*dv2-0.0625*B0[0]*favg[13]*dv2-0.0625*B0[3]*favg[12]*dv2-0.03608439182435161*B0[1]*favg[6]*dv2-0.03227486121839514*favg[2]*B0[5]*dv2-0.0625*B0[2]*favg[5]*dv2-0.03608439182435161*favg[1]*B0[3]*dv2-0.03608439182435161*favg[0]*B0[2]*dv2-0.03608439182435161*B0[0]*favg[2]*dv2+0.03608439182435161*B1[2]*favg[11]*dv1-0.25*fjump[10]*amax+0.2165063509461096*abar2[2]*favg[15]+0.1118033988749895*abar2[5]*favg[10]+0.125*abar2[0]*favg[10]+0.125*abar2[3]*favg[9]+0.125*abar2[2]*favg[4]; 
  incr[11] = (-0.0625*B0[0]*favg[14]*dv2)-0.03608439182435161*B0[2]*favg[8]*dv2-0.03608439182435161*B0[1]*favg[7]*dv2-0.03608439182435161*B0[0]*favg[3]*dv2+0.0625*B1[0]*favg[15]*dv1+0.03608439182435161*B1[2]*favg[10]*dv1+0.03608439182435161*B1[1]*favg[9]*dv1+0.03608439182435161*B1[0]*favg[4]*dv1-0.25*fjump[11]*amax+0.125*abar2[0]*favg[11]; 
  incr[12] = 0.1082531754730548*B0[1]*favg[15]*dv2+0.0625*B0[3]*favg[10]*dv2+0.05590169943749475*B0[4]*favg[9]*dv2+0.0625*B0[0]*favg[9]*dv2+0.0625*B0[1]*favg[4]*dv2-0.1082531754730548*B1[1]*favg[14]*dv1-0.0625*B1[3]*favg[8]*dv1-0.05590169943749475*B1[4]*favg[7]*dv1-0.0625*B1[0]*favg[7]*dv1-0.0625*B1[1]*favg[3]*dv1+0.75*fjump[12]*amax+0.4330127018922193*fjump[1]*amax-0.4841229182759272*abar2[1]*favg[20]-0.1936491673103709*abar2[1]*favg[16]-0.375*abar2[3]*favg[13]-0.3354101966249685*abar2[4]*favg[12]-0.375*abar2[0]*favg[12]-0.2165063509461096*abar2[2]*favg[6]-0.375*abar2[1]*favg[5]-0.1936491673103709*favg[1]*abar2[4]-0.2165063509461096*favg[2]*abar2[3]-0.2165063509461096*abar2[0]*favg[1]-0.2165063509461096*favg[0]*abar2[1]; 
  incr[13] = 0.1082531754730548*B0[2]*favg[15]*dv2+0.05590169943749475*B0[5]*favg[10]*dv2+0.0625*B0[0]*favg[10]*dv2+0.0625*B0[3]*favg[9]*dv2+0.0625*B0[2]*favg[4]*dv2-0.1082531754730548*B1[2]*favg[14]*dv1-0.05590169943749475*B1[5]*favg[8]*dv1-0.0625*B1[0]*favg[8]*dv1-0.0625*B1[3]*favg[7]*dv1-0.0625*B1[2]*favg[3]*dv1+0.75*fjump[13]*amax+0.4330127018922193*fjump[2]*amax-0.4841229182759272*abar2[2]*favg[20]-0.1936491673103709*abar2[2]*favg[17]-0.3354101966249685*abar2[5]*favg[13]-0.375*abar2[0]*favg[13]-0.375*abar2[3]*favg[12]-0.2165063509461096*abar2[1]*favg[6]-0.375*abar2[2]*favg[5]-0.1936491673103709*favg[2]*abar2[5]-0.2165063509461096*favg[1]*abar2[3]-0.2165063509461096*abar2[0]*favg[2]-0.2165063509461096*favg[0]*abar2[2]; 
  incr[14] = 0.0625*B0[0]*favg[11]*dv2-0.1397542485937369*B1[0]*favg[20]*dv1-0.05590169943749475*B1[0]*favg[18]*dv1-0.0625*B1[5]*favg[17]*dv1-0.0625*B1[4]*favg[16]*dv1-0.1082531754730548*B1[2]*favg[13]*dv1-0.1082531754730548*B1[1]*favg[12]*dv1-0.0625*B1[3]*favg[6]*dv1-0.1082531754730548*B1[0]*favg[5]*dv1-0.0625*favg[2]*B1[2]*dv1-0.0625*favg[1]*B1[1]*dv1-0.0625*favg[0]*B1[0]*dv1+0.75*fjump[14]*amax+0.4330127018922193*fjump[3]*amax-0.375*abar2[0]*favg[14]-0.2165063509461096*abar2[2]*favg[8]-0.2165063509461096*abar2[1]*favg[7]-0.2165063509461096*abar2[0]*favg[3]; 
  incr[15] = 0.1397542485937369*B0[0]*favg[20]*dv2+0.05590169943749475*B0[0]*favg[19]*dv2+0.0625*B0[5]*favg[17]*dv2+0.0625*B0[4]*favg[16]*dv2+0.1082531754730548*B0[2]*favg[13]*dv2+0.1082531754730548*B0[1]*favg[12]*dv2+0.0625*B0[3]*favg[6]*dv2+0.1082531754730548*B0[0]*favg[5]*dv2+0.0625*favg[2]*B0[2]*dv2+0.0625*favg[1]*B0[1]*dv2+0.0625*favg[0]*B0[0]*dv2-0.0625*B1[0]*favg[11]*dv1+0.75*fjump[15]*amax+0.4330127018922193*fjump[4]*amax-0.375*abar2[0]*favg[15]-0.2165063509461096*abar2[2]*favg[10]-0.2165063509461096*abar2[1]*favg[9]-0.2165063509461096*abar2[0]*favg[4]; 
  incr[16] = (-0.0625*B0[4]*favg[15]*dv2)-0.03227486121839514*B0[1]*favg[9]*dv2-0.03608439182435162*favg[4]*B0[4]*dv2+0.0625*B1[4]*favg[14]*dv1+0.03227486121839514*B1[1]*favg[7]*dv1+0.03608439182435162*favg[3]*B1[4]*dv1-0.25*fjump[16]*amax+0.2795084971874738*abar2[4]*favg[20]+0.0798595706249925*abar2[4]*favg[16]+0.125*abar2[0]*favg[16]+0.1936491673103709*abar2[1]*favg[12]+0.1118033988749895*abar2[3]*favg[6]+0.2165063509461096*abar2[4]*favg[5]+0.125*favg[0]*abar2[4]+0.1118033988749895*abar2[1]*favg[1]; 
  incr[17] = (-0.0625*B0[5]*favg[15]*dv2)-0.03227486121839514*B0[2]*favg[10]*dv2-0.03608439182435162*favg[4]*B0[5]*dv2+0.0625*B1[5]*favg[14]*dv1+0.03227486121839514*B1[2]*favg[8]*dv1+0.03608439182435162*favg[3]*B1[5]*dv1-0.25*fjump[17]*amax+0.2795084971874738*abar2[5]*favg[20]+0.0798595706249925*abar2[5]*favg[17]+0.125*abar2[0]*favg[17]+0.1936491673103709*abar2[2]*favg[13]+0.1118033988749895*abar2[3]*favg[6]+0.2165063509461096*abar2[5]*favg[5]+0.125*favg[0]*abar2[5]+0.1118033988749895*abar2[2]*favg[2]; 
  incr[18] = 0.05590169943749475*B1[0]*favg[14]*dv1+0.03227486121839514*B1[2]*favg[8]*dv1+0.03227486121839514*B1[1]*favg[7]*dv1+0.03227486121839514*B1[0]*favg[3]*dv1-0.25*fjump[18]*amax+0.125*abar2[0]*favg[18]; 
  incr[19] = (-0.05590169943749475*B0[0]*favg[15]*dv2)-0.03227486121839514*B0[2]*favg[10]*dv2-0.03227486121839514*B0[1]*favg[9]*dv2-0.03227486121839514*B0[0]*favg[4]*dv2-0.25*fjump[19]*amax+0.125*abar2[0]*favg[19]; 
  incr[20] = (-0.1397542485937369*B0[0]*favg[15]*dv2)-0.08068715304598785*B0[2]*favg[10]*dv2-0.08068715304598785*B0[1]*favg[9]*dv2-0.08068715304598785*B0[0]*favg[4]*dv2+0.1397542485937369*B1[0]*favg[14]*dv1+0.08068715304598785*B1[2]*favg[8]*dv1+0.08068715304598785*B1[1]*favg[7]*dv1+0.08068715304598785*B1[0]*favg[3]*dv1-1.25*fjump[20]*amax-0.9682458365518543*fjump[5]*amax-0.5590169943749475*fjump[0]*amax+0.625*abar2[0]*favg[20]+0.2795084971874737*abar2[5]*favg[17]+0.2795084971874737*abar2[4]*favg[16]+0.484122918275927*abar2[2]*favg[13]+0.484122918275927*abar2[1]*favg[12]+0.2795084971874737*abar2[3]*favg[6]+0.484122918275927*abar2[0]*favg[5]+0.2795084971874737*abar2[2]*favg[2]+0.2795084971874737*abar2[1]*favg[1]+0.2795084971874737*abar2[0]*favg[0]; 

  outr[0] += incr[0]*dv12; 
  outr[1] += incr[1]*dv12; 
  outr[2] += incr[2]*dv12; 
  outr[3] += incr[3]*dv12; 
  outr[4] += incr[4]*dv12; 
  outr[5] += incr[5]*dv12; 
  outr[6] += incr[6]*dv12; 
  outr[7] += incr[7]*dv12; 
  outr[8] += incr[8]*dv12; 
  outr[9] += incr[9]*dv12; 
  outr[10] += incr[10]*dv12; 
  outr[11] += incr[11]*dv12; 
  outr[12] += incr[12]*dv12; 
  outr[13] += incr[13]*dv12; 
  outr[14] += incr[14]*dv12; 
  outr[15] += incr[15]*dv12; 
  outr[16] += incr[16]*dv12; 
  outr[17] += incr[17]*dv12; 
  outr[18] += incr[18]*dv12; 
  outr[19] += incr[19]*dv12; 
  outr[20] += incr[20]*dv12; 

  outl[0] += -1.0*incr[0]*dv12; 
  outl[1] += -1.0*incr[1]*dv12; 
  outl[2] += -1.0*incr[2]*dv12; 
  outl[3] += -1.0*incr[3]*dv12; 
  outl[4] += -1.0*incr[4]*dv12; 
  outl[5] += incr[5]*dv12; 
  outl[6] += -1.0*incr[6]*dv12; 
  outl[7] += -1.0*incr[7]*dv12; 
  outl[8] += -1.0*incr[8]*dv12; 
  outl[9] += -1.0*incr[9]*dv12; 
  outl[10] += -1.0*incr[10]*dv12; 
  outl[11] += -1.0*incr[11]*dv12; 
  outl[12] += incr[12]*dv12; 
  outl[13] += incr[13]*dv12; 
  outl[14] += incr[14]*dv12; 
  outl[15] += incr[15]*dv12; 
  outl[16] += -1.0*incr[16]*dv12; 
  outl[17] += -1.0*incr[17]*dv12; 
  outl[18] += -1.0*incr[18]*dv12; 
  outl[19] += -1.0*incr[19]*dv12; 
  outl[20] += -1.0*incr[20]*dv12; 
return std::abs(amid); 
} 
