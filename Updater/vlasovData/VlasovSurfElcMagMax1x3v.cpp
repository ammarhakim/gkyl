#include <VlasovModDecl.h> 
double VlasovSurfElcMag1x3vMax_VX_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 

  double abar2[2]; 


  abar0[0] = E0[0]+wv2*B2[0]-wv3*B1[0]; 
  abar0[1] = E0[1]+wv2*B2[1]-wv3*B1[1]; 

  double incr[5]; 

  double favg[5]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  double fjump[5]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  const double amid = 0.7071067811865475*abar0[0]; 
  incr[0] = (-0.05103103630798287*B1[0]*favg[4]*dv3)+0.05103103630798287*B2[0]*favg[3]*dv2-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.05103103630798287*B1[1]*favg[4]*dv3)+0.05103103630798287*B2[1]*favg[3]*dv2-0.25*fjump[1]*amax+0.3061862178478973*abar0[1]*favg[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.08838834764831845*B1[0]*favg[4]*dv3-0.08838834764831845*B2[0]*favg[3]*dv2+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = 0.08838834764831845*B2[0]*favg[2]*dv2+0.05103103630798287*favg[1]*B2[1]*dv2+0.05103103630798287*favg[0]*B2[0]*dv2-0.25*fjump[3]*amax+0.1767766952966369*abar0[0]*favg[3]; 
  incr[4] = (-0.08838834764831845*B1[0]*favg[2]*dv3)-0.05103103630798287*favg[1]*B1[1]*dv3-0.05103103630798287*favg[0]*B1[0]*dv3-0.25*fjump[4]*amax+0.1767766952966369*abar0[0]*favg[4]; 

  outr[0] += incr[0]*dv10; 
  outr[1] += incr[1]*dv10; 
  outr[2] += incr[2]*dv10; 
  outr[3] += incr[3]*dv10; 
  outr[4] += incr[4]*dv10; 

  outl[0] += -1.0*incr[0]*dv10; 
  outl[1] += -1.0*incr[1]*dv10; 
  outl[2] += incr[2]*dv10; 
  outl[3] += -1.0*incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x3vMax_VX_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar0[0] = E0[0]+wv2*B2[0]-wv3*B1[0]; 
  abar0[1] = E0[1]+wv2*B2[1]-wv3*B1[1]; 
  abar0[2] = E0[2]+wv2*B2[2]-wv3*B1[2]; 

  double incr[15]; 

  double favg[15]; 

  favg[0] = 1*fr[0]+fl[0]; 
  favg[1] = 1*fr[1]+fl[1]; 
  favg[2] = -1*fr[2]+fl[2]; 
  favg[3] = 1*fr[3]+fl[3]; 
  favg[4] = 1*fr[4]+fl[4]; 
  favg[5] = -1*fr[5]+fl[5]; 
  favg[6] = 1*fr[6]+fl[6]; 
  favg[7] = -1*fr[7]+fl[7]; 
  favg[8] = 1*fr[8]+fl[8]; 
  favg[9] = -1*fr[9]+fl[9]; 
  favg[10] = 1*fr[10]+fl[10]; 
  favg[11] = 1*fr[11]+fl[11]; 
  favg[12] = 1*fr[12]+fl[12]; 
  favg[13] = 1*fr[13]+fl[13]; 
  favg[14] = 1*fr[14]+fl[14]; 
  double fjump[15]; 

  fjump[0] = 1*fr[0]-fl[0]; 
  fjump[1] = 1*fr[1]-fl[1]; 
  fjump[2] = -1*fr[2]-fl[2]; 
  fjump[3] = 1*fr[3]-fl[3]; 
  fjump[4] = 1*fr[4]-fl[4]; 
  fjump[5] = -1*fr[5]-fl[5]; 
  fjump[6] = 1*fr[6]-fl[6]; 
  fjump[7] = -1*fr[7]-fl[7]; 
  fjump[8] = 1*fr[8]-fl[8]; 
  fjump[9] = -1*fr[9]-fl[9]; 
  fjump[10] = 1*fr[10]-fl[10]; 
  fjump[11] = 1*fr[11]-fl[11]; 
  fjump[12] = 1*fr[12]-fl[12]; 
  fjump[13] = 1*fr[13]-fl[13]; 
  fjump[14] = 1*fr[14]-fl[14]; 
  const double amid = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr[0] = (-0.08838834764831845*B1[0]*favg[9]*dv3)-0.05103103630798287*B1[1]*favg[8]*dv3-0.05103103630798287*B1[0]*favg[4]*dv3+0.08838834764831845*B2[0]*favg[7]*dv2+0.05103103630798287*B2[1]*favg[6]*dv2+0.05103103630798287*B2[0]*favg[3]*dv2-0.5590169943749475*fjump[12]*amax-0.4330127018922193*fjump[2]*amax-0.25*fjump[0]*amax+0.3952847075210476*abar0[0]*favg[12]+0.1767766952966369*abar0[2]*favg[11]+0.3061862178478973*abar0[1]*favg[5]+0.3061862178478973*abar0[0]*favg[2]+0.1767766952966369*abar0[1]*favg[1]+0.1767766952966369*abar0[0]*favg[0]; 
  incr[1] = (-0.08838834764831845*B1[1]*favg[9]*dv3)-0.04564354645876385*B1[2]*favg[8]*dv3-0.05103103630798287*B1[0]*favg[8]*dv3-0.05103103630798287*B1[1]*favg[4]*dv3+0.08838834764831845*B2[1]*favg[7]*dv2+0.04564354645876385*B2[2]*favg[6]*dv2+0.05103103630798287*B2[0]*favg[6]*dv2+0.05103103630798287*B2[1]*favg[3]*dv2-0.4330127018922193*fjump[5]*amax-0.25*fjump[1]*amax+0.3952847075210476*abar0[1]*favg[12]+0.158113883008419*abar0[1]*favg[11]+0.2738612787525831*abar0[2]*favg[5]+0.3061862178478973*abar0[0]*favg[5]+0.3061862178478973*abar0[1]*favg[2]+0.158113883008419*favg[1]*abar0[2]+0.1767766952966369*abar0[0]*favg[1]+0.1767766952966369*favg[0]*abar0[1]; 
  incr[2] = 0.1530931089239486*B1[0]*favg[9]*dv3+0.08838834764831845*B1[1]*favg[8]*dv3+0.08838834764831845*B1[0]*favg[4]*dv3-0.1530931089239486*B2[0]*favg[7]*dv2-0.08838834764831845*B2[1]*favg[6]*dv2-0.08838834764831845*B2[0]*favg[3]*dv2+0.9682458365518543*fjump[12]*amax+0.75*fjump[2]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814578*abar0[0]*favg[12]-0.3061862178478973*abar0[2]*favg[11]-0.5303300858899107*abar0[1]*favg[5]-0.5303300858899107*abar0[0]*favg[2]-0.3061862178478973*abar0[1]*favg[1]-0.3061862178478973*abar0[0]*favg[0]; 
  incr[3] = (-0.05103103630798287*B1[0]*favg[10]*dv3)+0.04564354645876385*B2[0]*favg[13]*dv2+0.1141088661469096*B2[0]*favg[12]*dv2+0.05103103630798287*B2[2]*favg[11]*dv2+0.08838834764831845*B2[1]*favg[5]*dv2+0.08838834764831845*B2[0]*favg[2]*dv2+0.05103103630798287*favg[1]*B2[1]*dv2+0.05103103630798287*favg[0]*B2[0]*dv2-0.4330127018922193*fjump[7]*amax-0.25*fjump[3]*amax+0.3061862178478973*abar0[0]*favg[7]+0.1767766952966369*abar0[1]*favg[6]+0.1767766952966369*abar0[0]*favg[3]; 
  incr[4] = (-0.04564354645876385*B1[0]*favg[14]*dv3)-0.1141088661469096*B1[0]*favg[12]*dv3-0.05103103630798287*B1[2]*favg[11]*dv3-0.08838834764831845*B1[1]*favg[5]*dv3-0.08838834764831845*B1[0]*favg[2]*dv3-0.05103103630798287*favg[1]*B1[1]*dv3-0.05103103630798287*favg[0]*B1[0]*dv3+0.05103103630798287*B2[0]*favg[10]*dv2-0.4330127018922193*fjump[9]*amax-0.25*fjump[4]*amax+0.3061862178478973*abar0[0]*favg[9]+0.1767766952966369*abar0[1]*favg[8]+0.1767766952966369*abar0[0]*favg[4]; 
  incr[5] = 0.1530931089239486*B1[1]*favg[9]*dv3+0.0790569415042095*B1[2]*favg[8]*dv3+0.08838834764831845*B1[0]*favg[8]*dv3+0.08838834764831845*B1[1]*favg[4]*dv3-0.1530931089239486*B2[1]*favg[7]*dv2-0.0790569415042095*B2[2]*favg[6]*dv2-0.08838834764831845*B2[0]*favg[6]*dv2-0.08838834764831845*B2[1]*favg[3]*dv2+0.75*fjump[5]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814579*abar0[1]*favg[12]-0.2738612787525831*abar0[1]*favg[11]-0.4743416490252572*abar0[2]*favg[5]-0.5303300858899107*abar0[0]*favg[5]-0.5303300858899107*abar0[1]*favg[2]-0.2738612787525831*favg[1]*abar0[2]-0.3061862178478973*abar0[0]*favg[1]-0.3061862178478973*favg[0]*abar0[1]; 
  incr[6] = (-0.05103103630798287*B1[1]*favg[10]*dv3)+0.04564354645876385*B2[1]*favg[13]*dv2+0.1141088661469096*B2[1]*favg[12]*dv2+0.04564354645876385*B2[1]*favg[11]*dv2+0.0790569415042095*B2[2]*favg[5]*dv2+0.08838834764831845*B2[0]*favg[5]*dv2+0.04564354645876385*favg[1]*B2[2]*dv2+0.08838834764831845*B2[1]*favg[2]*dv2+0.05103103630798287*favg[0]*B2[1]*dv2+0.05103103630798287*B2[0]*favg[1]*dv2-0.25*fjump[6]*amax+0.3061862178478973*abar0[1]*favg[7]+0.158113883008419*abar0[2]*favg[6]+0.1767766952966369*abar0[0]*favg[6]+0.1767766952966369*abar0[1]*favg[3]; 
  incr[7] = 0.08838834764831845*B1[0]*favg[10]*dv3-0.0790569415042095*B2[0]*favg[13]*dv2-0.1976423537605238*B2[0]*favg[12]*dv2-0.08838834764831845*B2[2]*favg[11]*dv2-0.1530931089239486*B2[1]*favg[5]*dv2-0.1530931089239486*B2[0]*favg[2]*dv2-0.08838834764831845*favg[1]*B2[1]*dv2-0.08838834764831845*favg[0]*B2[0]*dv2+0.75*fjump[7]*amax+0.4330127018922193*fjump[3]*amax-0.5303300858899107*abar0[0]*favg[7]-0.3061862178478973*abar0[1]*favg[6]-0.3061862178478973*abar0[0]*favg[3]; 
  incr[8] = (-0.04564354645876385*B1[1]*favg[14]*dv3)-0.1141088661469096*B1[1]*favg[12]*dv3-0.04564354645876385*B1[1]*favg[11]*dv3-0.0790569415042095*B1[2]*favg[5]*dv3-0.08838834764831845*B1[0]*favg[5]*dv3-0.04564354645876385*favg[1]*B1[2]*dv3-0.08838834764831845*B1[1]*favg[2]*dv3-0.05103103630798287*favg[0]*B1[1]*dv3-0.05103103630798287*B1[0]*favg[1]*dv3+0.05103103630798287*B2[1]*favg[10]*dv2-0.25*fjump[8]*amax+0.3061862178478973*abar0[1]*favg[9]+0.158113883008419*abar0[2]*favg[8]+0.1767766952966369*abar0[0]*favg[8]+0.1767766952966369*abar0[1]*favg[4]; 
  incr[9] = 0.0790569415042095*B1[0]*favg[14]*dv3+0.1976423537605238*B1[0]*favg[12]*dv3+0.08838834764831845*B1[2]*favg[11]*dv3+0.1530931089239486*B1[1]*favg[5]*dv3+0.1530931089239486*B1[0]*favg[2]*dv3+0.08838834764831845*favg[1]*B1[1]*dv3+0.08838834764831845*favg[0]*B1[0]*dv3-0.08838834764831845*B2[0]*favg[10]*dv2+0.75*fjump[9]*amax+0.4330127018922193*fjump[4]*amax-0.5303300858899107*abar0[0]*favg[9]-0.3061862178478973*abar0[1]*favg[8]-0.3061862178478973*abar0[0]*favg[4]; 
  incr[10] = (-0.08838834764831845*B1[0]*favg[7]*dv3)-0.05103103630798287*B1[1]*favg[6]*dv3-0.05103103630798287*B1[0]*favg[3]*dv3+0.08838834764831845*B2[0]*favg[9]*dv2+0.05103103630798287*B2[1]*favg[8]*dv2+0.05103103630798287*B2[0]*favg[4]*dv2-0.25*fjump[10]*amax+0.1767766952966369*abar0[0]*favg[10]; 
  incr[11] = (-0.08838834764831845*B1[2]*favg[9]*dv3)-0.04564354645876386*B1[1]*favg[8]*dv3-0.05103103630798288*B1[2]*favg[4]*dv3+0.08838834764831845*B2[2]*favg[7]*dv2+0.04564354645876386*B2[1]*favg[6]*dv2+0.05103103630798288*B2[2]*favg[3]*dv2-0.25*fjump[11]*amax+0.3952847075210476*abar0[2]*favg[12]+0.1129384878631565*abar0[2]*favg[11]+0.1767766952966369*abar0[0]*favg[11]+0.2738612787525832*abar0[1]*favg[5]+0.3061862178478973*abar0[2]*favg[2]+0.1767766952966369*favg[0]*abar0[2]+0.158113883008419*abar0[1]*favg[1]; 
  incr[12] = (-0.1976423537605237*B1[0]*favg[9]*dv3)-0.1141088661469096*B1[1]*favg[8]*dv3-0.1141088661469096*B1[0]*favg[4]*dv3+0.1976423537605237*B2[0]*favg[7]*dv2+0.1141088661469096*B2[1]*favg[6]*dv2+0.1141088661469096*B2[0]*favg[3]*dv2-1.25*fjump[12]*amax-0.9682458365518543*fjump[2]*amax-0.5590169943749475*fjump[0]*amax+0.8838834764831847*abar0[0]*favg[12]+0.3952847075210476*abar0[2]*favg[11]+0.6846531968814578*abar0[1]*favg[5]+0.6846531968814578*abar0[0]*favg[2]+0.3952847075210476*abar0[1]*favg[1]+0.3952847075210476*abar0[0]*favg[0]; 
  incr[13] = 0.0790569415042095*B2[0]*favg[7]*dv2+0.04564354645876384*B2[1]*favg[6]*dv2+0.04564354645876384*B2[0]*favg[3]*dv2-0.25*fjump[13]*amax+0.1767766952966369*abar0[0]*favg[13]; 
  incr[14] = (-0.0790569415042095*B1[0]*favg[9]*dv3)-0.04564354645876384*B1[1]*favg[8]*dv3-0.04564354645876384*B1[0]*favg[4]*dv3-0.25*fjump[14]*amax+0.1767766952966369*abar0[0]*favg[14]; 

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
  outl[2] += incr[2]*dv10; 
  outl[3] += -1.0*incr[3]*dv10; 
  outl[4] += -1.0*incr[4]*dv10; 
  outl[5] += incr[5]*dv10; 
  outl[6] += -1.0*incr[6]*dv10; 
  outl[7] += incr[7]*dv10; 
  outl[8] += -1.0*incr[8]*dv10; 
  outl[9] += incr[9]*dv10; 
  outl[10] += -1.0*incr[10]*dv10; 
  outl[11] += -1.0*incr[11]*dv10; 
  outl[12] += -1.0*incr[12]*dv10; 
  outl[13] += -1.0*incr[13]*dv10; 
  outl[14] += -1.0*incr[14]*dv10; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x3vMax_VY_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[2]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 

  double abar2[2]; 


  abar1[0] = E1[0]+wv3*B0[0]-wv1*B2[0]; 
  abar1[1] = E1[1]+wv3*B0[1]-wv1*B2[1]; 

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
  const double amid = 0.7071067811865475*abar1[0]; 
  incr[0] = 0.05103103630798287*B0[0]*favg[4]*dv3-0.05103103630798287*B2[0]*favg[2]*dv1-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.3061862178478973*abar1[0]*favg[3]+0.1767766952966369*abar1[1]*favg[1]+0.1767766952966369*abar1[0]*favg[0]; 
  incr[1] = 0.05103103630798287*B0[1]*favg[4]*dv3-0.05103103630798287*B2[1]*favg[2]*dv1-0.25*fjump[1]*amax+0.3061862178478973*abar1[1]*favg[3]+0.1767766952966369*abar1[0]*favg[1]+0.1767766952966369*favg[0]*abar1[1]; 
  incr[2] = (-0.08838834764831845*B2[0]*favg[3]*dv1)-0.05103103630798287*favg[1]*B2[1]*dv1-0.05103103630798287*favg[0]*B2[0]*dv1-0.25*fjump[2]*amax+0.1767766952966369*abar1[0]*favg[2]; 
  incr[3] = (-0.08838834764831845*B0[0]*favg[4]*dv3)+0.08838834764831845*B2[0]*favg[2]*dv1+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899107*abar1[0]*favg[3]-0.3061862178478973*abar1[1]*favg[1]-0.3061862178478973*abar1[0]*favg[0]; 
  incr[4] = 0.08838834764831845*B0[0]*favg[3]*dv3+0.05103103630798287*favg[1]*B0[1]*dv3+0.05103103630798287*favg[0]*B0[0]*dv3-0.25*fjump[4]*amax+0.1767766952966369*abar1[0]*favg[4]; 

  outr[0] += incr[0]*dv11; 
  outr[1] += incr[1]*dv11; 
  outr[2] += incr[2]*dv11; 
  outr[3] += incr[3]*dv11; 
  outr[4] += incr[4]*dv11; 

  outl[0] += -1.0*incr[0]*dv11; 
  outl[1] += -1.0*incr[1]*dv11; 
  outl[2] += -1.0*incr[2]*dv11; 
  outl[3] += incr[3]*dv11; 
  outl[4] += -1.0*incr[4]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x3vMax_VY_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[3]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar1[0] = E1[0]+wv3*B0[0]-wv1*B2[0]; 
  abar1[1] = E1[1]+wv3*B0[1]-wv1*B2[1]; 
  abar1[2] = E1[2]+wv3*B0[2]-wv1*B2[2]; 

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
  const double amid = 0.7071067811865475*abar1[0]-0.7905694150420947*abar1[2]; 
  incr[0] = 0.08838834764831845*B0[0]*favg[10]*dv3+0.05103103630798287*B0[1]*favg[8]*dv3+0.05103103630798287*B0[0]*favg[4]*dv3-0.08838834764831845*B2[0]*favg[7]*dv1-0.05103103630798287*B2[1]*favg[5]*dv1-0.05103103630798287*B2[0]*favg[2]*dv1-0.5590169943749475*fjump[13]*amax-0.4330127018922193*fjump[3]*amax-0.25*fjump[0]*amax+0.3952847075210476*abar1[0]*favg[13]+0.1767766952966369*abar1[2]*favg[11]+0.3061862178478973*abar1[1]*favg[6]+0.3061862178478973*abar1[0]*favg[3]+0.1767766952966369*abar1[1]*favg[1]+0.1767766952966369*abar1[0]*favg[0]; 
  incr[1] = 0.08838834764831845*B0[1]*favg[10]*dv3+0.04564354645876385*B0[2]*favg[8]*dv3+0.05103103630798287*B0[0]*favg[8]*dv3+0.05103103630798287*B0[1]*favg[4]*dv3-0.08838834764831845*B2[1]*favg[7]*dv1-0.04564354645876385*B2[2]*favg[5]*dv1-0.05103103630798287*B2[0]*favg[5]*dv1-0.05103103630798287*B2[1]*favg[2]*dv1-0.4330127018922193*fjump[6]*amax-0.25*fjump[1]*amax+0.3952847075210476*abar1[1]*favg[13]+0.158113883008419*abar1[1]*favg[11]+0.2738612787525831*abar1[2]*favg[6]+0.3061862178478973*abar1[0]*favg[6]+0.3061862178478973*abar1[1]*favg[3]+0.158113883008419*favg[1]*abar1[2]+0.1767766952966369*abar1[0]*favg[1]+0.1767766952966369*favg[0]*abar1[1]; 
  incr[2] = 0.05103103630798287*B0[0]*favg[9]*dv3-0.1141088661469096*B2[0]*favg[13]*dv1-0.04564354645876385*B2[0]*favg[12]*dv1-0.05103103630798287*B2[2]*favg[11]*dv1-0.08838834764831845*B2[1]*favg[6]*dv1-0.08838834764831845*B2[0]*favg[3]*dv1-0.05103103630798287*favg[1]*B2[1]*dv1-0.05103103630798287*favg[0]*B2[0]*dv1-0.4330127018922193*fjump[7]*amax-0.25*fjump[2]*amax+0.3061862178478973*abar1[0]*favg[7]+0.1767766952966369*abar1[1]*favg[5]+0.1767766952966369*abar1[0]*favg[2]; 
  incr[3] = (-0.1530931089239486*B0[0]*favg[10]*dv3)-0.08838834764831845*B0[1]*favg[8]*dv3-0.08838834764831845*B0[0]*favg[4]*dv3+0.1530931089239486*B2[0]*favg[7]*dv1+0.08838834764831845*B2[1]*favg[5]*dv1+0.08838834764831845*B2[0]*favg[2]*dv1+0.9682458365518543*fjump[13]*amax+0.75*fjump[3]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814578*abar1[0]*favg[13]-0.3061862178478973*abar1[2]*favg[11]-0.5303300858899107*abar1[1]*favg[6]-0.5303300858899107*abar1[0]*favg[3]-0.3061862178478973*abar1[1]*favg[1]-0.3061862178478973*abar1[0]*favg[0]; 
  incr[4] = 0.04564354645876385*B0[0]*favg[14]*dv3+0.1141088661469096*B0[0]*favg[13]*dv3+0.05103103630798287*B0[2]*favg[11]*dv3+0.08838834764831845*B0[1]*favg[6]*dv3+0.08838834764831845*B0[0]*favg[3]*dv3+0.05103103630798287*favg[1]*B0[1]*dv3+0.05103103630798287*favg[0]*B0[0]*dv3-0.05103103630798287*B2[0]*favg[9]*dv1-0.4330127018922193*fjump[10]*amax-0.25*fjump[4]*amax+0.3061862178478973*abar1[0]*favg[10]+0.1767766952966369*abar1[1]*favg[8]+0.1767766952966369*abar1[0]*favg[4]; 
  incr[5] = 0.05103103630798287*B0[1]*favg[9]*dv3-0.1141088661469096*B2[1]*favg[13]*dv1-0.04564354645876385*B2[1]*favg[12]*dv1-0.04564354645876385*B2[1]*favg[11]*dv1-0.0790569415042095*B2[2]*favg[6]*dv1-0.08838834764831845*B2[0]*favg[6]*dv1-0.08838834764831845*B2[1]*favg[3]*dv1-0.04564354645876385*favg[1]*B2[2]*dv1-0.05103103630798287*favg[0]*B2[1]*dv1-0.05103103630798287*B2[0]*favg[1]*dv1-0.25*fjump[5]*amax+0.3061862178478973*abar1[1]*favg[7]+0.158113883008419*abar1[2]*favg[5]+0.1767766952966369*abar1[0]*favg[5]+0.1767766952966369*abar1[1]*favg[2]; 
  incr[6] = (-0.1530931089239486*B0[1]*favg[10]*dv3)-0.0790569415042095*B0[2]*favg[8]*dv3-0.08838834764831845*B0[0]*favg[8]*dv3-0.08838834764831845*B0[1]*favg[4]*dv3+0.1530931089239486*B2[1]*favg[7]*dv1+0.0790569415042095*B2[2]*favg[5]*dv1+0.08838834764831845*B2[0]*favg[5]*dv1+0.08838834764831845*B2[1]*favg[2]*dv1+0.75*fjump[6]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814579*abar1[1]*favg[13]-0.2738612787525831*abar1[1]*favg[11]-0.4743416490252572*abar1[2]*favg[6]-0.5303300858899107*abar1[0]*favg[6]-0.5303300858899107*abar1[1]*favg[3]-0.2738612787525831*favg[1]*abar1[2]-0.3061862178478973*abar1[0]*favg[1]-0.3061862178478973*favg[0]*abar1[1]; 
  incr[7] = (-0.08838834764831845*B0[0]*favg[9]*dv3)+0.1976423537605238*B2[0]*favg[13]*dv1+0.0790569415042095*B2[0]*favg[12]*dv1+0.08838834764831845*B2[2]*favg[11]*dv1+0.1530931089239486*B2[1]*favg[6]*dv1+0.1530931089239486*B2[0]*favg[3]*dv1+0.08838834764831845*favg[1]*B2[1]*dv1+0.08838834764831845*favg[0]*B2[0]*dv1+0.75*fjump[7]*amax+0.4330127018922193*fjump[2]*amax-0.5303300858899107*abar1[0]*favg[7]-0.3061862178478973*abar1[1]*favg[5]-0.3061862178478973*abar1[0]*favg[2]; 
  incr[8] = 0.04564354645876385*B0[1]*favg[14]*dv3+0.1141088661469096*B0[1]*favg[13]*dv3+0.04564354645876385*B0[1]*favg[11]*dv3+0.0790569415042095*B0[2]*favg[6]*dv3+0.08838834764831845*B0[0]*favg[6]*dv3+0.08838834764831845*B0[1]*favg[3]*dv3+0.04564354645876385*favg[1]*B0[2]*dv3+0.05103103630798287*favg[0]*B0[1]*dv3+0.05103103630798287*B0[0]*favg[1]*dv3-0.05103103630798287*B2[1]*favg[9]*dv1-0.25*fjump[8]*amax+0.3061862178478973*abar1[1]*favg[10]+0.158113883008419*abar1[2]*favg[8]+0.1767766952966369*abar1[0]*favg[8]+0.1767766952966369*abar1[1]*favg[4]; 
  incr[9] = 0.08838834764831845*B0[0]*favg[7]*dv3+0.05103103630798287*B0[1]*favg[5]*dv3+0.05103103630798287*B0[0]*favg[2]*dv3-0.08838834764831845*B2[0]*favg[10]*dv1-0.05103103630798287*B2[1]*favg[8]*dv1-0.05103103630798287*B2[0]*favg[4]*dv1-0.25*fjump[9]*amax+0.1767766952966369*abar1[0]*favg[9]; 
  incr[10] = (-0.0790569415042095*B0[0]*favg[14]*dv3)-0.1976423537605238*B0[0]*favg[13]*dv3-0.08838834764831845*B0[2]*favg[11]*dv3-0.1530931089239486*B0[1]*favg[6]*dv3-0.1530931089239486*B0[0]*favg[3]*dv3-0.08838834764831845*favg[1]*B0[1]*dv3-0.08838834764831845*favg[0]*B0[0]*dv3+0.08838834764831845*B2[0]*favg[9]*dv1+0.75*fjump[10]*amax+0.4330127018922193*fjump[4]*amax-0.5303300858899107*abar1[0]*favg[10]-0.3061862178478973*abar1[1]*favg[8]-0.3061862178478973*abar1[0]*favg[4]; 
  incr[11] = 0.08838834764831845*B0[2]*favg[10]*dv3+0.04564354645876386*B0[1]*favg[8]*dv3+0.05103103630798288*B0[2]*favg[4]*dv3-0.08838834764831845*B2[2]*favg[7]*dv1-0.04564354645876386*B2[1]*favg[5]*dv1-0.05103103630798288*favg[2]*B2[2]*dv1-0.25*fjump[11]*amax+0.3952847075210476*abar1[2]*favg[13]+0.1129384878631565*abar1[2]*favg[11]+0.1767766952966369*abar1[0]*favg[11]+0.2738612787525832*abar1[1]*favg[6]+0.3061862178478973*abar1[2]*favg[3]+0.1767766952966369*favg[0]*abar1[2]+0.158113883008419*abar1[1]*favg[1]; 
  incr[12] = (-0.0790569415042095*B2[0]*favg[7]*dv1)-0.04564354645876384*B2[1]*favg[5]*dv1-0.04564354645876384*B2[0]*favg[2]*dv1-0.25*fjump[12]*amax+0.1767766952966369*abar1[0]*favg[12]; 
  incr[13] = 0.1976423537605237*B0[0]*favg[10]*dv3+0.1141088661469096*B0[1]*favg[8]*dv3+0.1141088661469096*B0[0]*favg[4]*dv3-0.1976423537605237*B2[0]*favg[7]*dv1-0.1141088661469096*B2[1]*favg[5]*dv1-0.1141088661469096*B2[0]*favg[2]*dv1-1.25*fjump[13]*amax-0.9682458365518543*fjump[3]*amax-0.5590169943749475*fjump[0]*amax+0.8838834764831847*abar1[0]*favg[13]+0.3952847075210476*abar1[2]*favg[11]+0.6846531968814578*abar1[1]*favg[6]+0.6846531968814578*abar1[0]*favg[3]+0.3952847075210476*abar1[1]*favg[1]+0.3952847075210476*abar1[0]*favg[0]; 
  incr[14] = 0.0790569415042095*B0[0]*favg[10]*dv3+0.04564354645876384*B0[1]*favg[8]*dv3+0.04564354645876384*B0[0]*favg[4]*dv3-0.25*fjump[14]*amax+0.1767766952966369*abar1[0]*favg[14]; 

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
  outl[3] += incr[3]*dv11; 
  outl[4] += -1.0*incr[4]*dv11; 
  outl[5] += -1.0*incr[5]*dv11; 
  outl[6] += incr[6]*dv11; 
  outl[7] += incr[7]*dv11; 
  outl[8] += -1.0*incr[8]*dv11; 
  outl[9] += -1.0*incr[9]*dv11; 
  outl[10] += incr[10]*dv11; 
  outl[11] += -1.0*incr[11]*dv11; 
  outl[12] += -1.0*incr[12]*dv11; 
  outl[13] += -1.0*incr[13]*dv11; 
  outl[14] += -1.0*incr[14]*dv11; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x3vMax_VZ_P1(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12 = 2/dxv[3]; 
  const double *E2 = &EM[4]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 

  double abar2[2]; 


  abar2[0] = E2[0]+wv1*B1[0]-wv2*B0[0]; 
  abar2[1] = E2[1]+wv1*B1[1]-wv2*B0[1]; 

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
  const double amid = 0.7071067811865475*abar2[0]; 
  incr[0] = (-0.05103103630798287*B0[0]*favg[3]*dv2)+0.05103103630798287*B1[0]*favg[2]*dv1-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.3061862178478973*abar2[0]*favg[4]+0.1767766952966369*abar2[1]*favg[1]+0.1767766952966369*abar2[0]*favg[0]; 
  incr[1] = (-0.05103103630798287*B0[1]*favg[3]*dv2)+0.05103103630798287*B1[1]*favg[2]*dv1-0.25*fjump[1]*amax+0.3061862178478973*abar2[1]*favg[4]+0.1767766952966369*abar2[0]*favg[1]+0.1767766952966369*favg[0]*abar2[1]; 
  incr[2] = 0.08838834764831845*B1[0]*favg[4]*dv1+0.05103103630798287*favg[1]*B1[1]*dv1+0.05103103630798287*favg[0]*B1[0]*dv1-0.25*fjump[2]*amax+0.1767766952966369*abar2[0]*favg[2]; 
  incr[3] = (-0.08838834764831845*B0[0]*favg[4]*dv2)-0.05103103630798287*favg[1]*B0[1]*dv2-0.05103103630798287*favg[0]*B0[0]*dv2-0.25*fjump[3]*amax+0.1767766952966369*abar2[0]*favg[3]; 
  incr[4] = 0.08838834764831845*B0[0]*favg[3]*dv2-0.08838834764831845*B1[0]*favg[2]*dv1+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.5303300858899107*abar2[0]*favg[4]-0.3061862178478973*abar2[1]*favg[1]-0.3061862178478973*abar2[0]*favg[0]; 

  outr[0] += incr[0]*dv12; 
  outr[1] += incr[1]*dv12; 
  outr[2] += incr[2]*dv12; 
  outr[3] += incr[3]*dv12; 
  outr[4] += incr[4]*dv12; 

  outl[0] += -1.0*incr[0]*dv12; 
  outl[1] += -1.0*incr[1]*dv12; 
  outl[2] += -1.0*incr[2]*dv12; 
  outl[3] += -1.0*incr[3]*dv12; 
  outl[4] += incr[4]*dv12; 
return std::abs(amid); 
} 
double VlasovSurfElcMag1x3vMax_VZ_P2(const double *w, const double *dxv, const double amax, const double *EM, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w: Cell-center coordinates. dxv[NDIM]: Cell spacing. amax: amax in global lax flux. E: EM field. fl/fr: Distribution function in left/right cells 
// outl/outr: output distribution function in left/right cells 
// returns abs(amid) for use in determining amax in cfl and global lax flux 
  double dv12 = 2/dxv[3]; 
  const double *E2 = &EM[6]; 

  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 

  double abar2[3]; 


  abar2[0] = E2[0]+wv1*B1[0]-wv2*B0[0]; 
  abar2[1] = E2[1]+wv1*B1[1]-wv2*B0[1]; 
  abar2[2] = E2[2]+wv1*B1[2]-wv2*B0[2]; 

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
  const double amid = 0.7071067811865475*abar2[0]-0.7905694150420947*abar2[2]; 
  incr[0] = (-0.08838834764831845*B0[0]*favg[10]*dv2)-0.05103103630798287*B0[1]*favg[6]*dv2-0.05103103630798287*B0[0]*favg[3]*dv2+0.08838834764831845*B1[0]*favg[9]*dv1+0.05103103630798287*B1[1]*favg[5]*dv1+0.05103103630798287*B1[0]*favg[2]*dv1-0.5590169943749475*fjump[14]*amax-0.4330127018922193*fjump[4]*amax-0.25*fjump[0]*amax+0.3952847075210476*abar2[0]*favg[14]+0.1767766952966369*abar2[2]*favg[11]+0.3061862178478973*abar2[1]*favg[8]+0.3061862178478973*abar2[0]*favg[4]+0.1767766952966369*abar2[1]*favg[1]+0.1767766952966369*abar2[0]*favg[0]; 
  incr[1] = (-0.08838834764831845*B0[1]*favg[10]*dv2)-0.04564354645876385*B0[2]*favg[6]*dv2-0.05103103630798287*B0[0]*favg[6]*dv2-0.05103103630798287*B0[1]*favg[3]*dv2+0.08838834764831845*B1[1]*favg[9]*dv1+0.04564354645876385*B1[2]*favg[5]*dv1+0.05103103630798287*B1[0]*favg[5]*dv1+0.05103103630798287*B1[1]*favg[2]*dv1-0.4330127018922193*fjump[8]*amax-0.25*fjump[1]*amax+0.3952847075210476*abar2[1]*favg[14]+0.158113883008419*abar2[1]*favg[11]+0.2738612787525831*abar2[2]*favg[8]+0.3061862178478973*abar2[0]*favg[8]+0.3061862178478973*abar2[1]*favg[4]+0.158113883008419*favg[1]*abar2[2]+0.1767766952966369*abar2[0]*favg[1]+0.1767766952966369*favg[0]*abar2[1]; 
  incr[2] = (-0.05103103630798287*B0[0]*favg[7]*dv2)+0.1141088661469096*B1[0]*favg[14]*dv1+0.04564354645876385*B1[0]*favg[12]*dv1+0.05103103630798287*B1[2]*favg[11]*dv1+0.08838834764831845*B1[1]*favg[8]*dv1+0.08838834764831845*B1[0]*favg[4]*dv1+0.05103103630798287*favg[1]*B1[1]*dv1+0.05103103630798287*favg[0]*B1[0]*dv1-0.4330127018922193*fjump[9]*amax-0.25*fjump[2]*amax+0.3061862178478973*abar2[0]*favg[9]+0.1767766952966369*abar2[1]*favg[5]+0.1767766952966369*abar2[0]*favg[2]; 
  incr[3] = (-0.1141088661469096*B0[0]*favg[14]*dv2)-0.04564354645876385*B0[0]*favg[13]*dv2-0.05103103630798287*B0[2]*favg[11]*dv2-0.08838834764831845*B0[1]*favg[8]*dv2-0.08838834764831845*B0[0]*favg[4]*dv2-0.05103103630798287*favg[1]*B0[1]*dv2-0.05103103630798287*favg[0]*B0[0]*dv2+0.05103103630798287*B1[0]*favg[7]*dv1-0.4330127018922193*fjump[10]*amax-0.25*fjump[3]*amax+0.3061862178478973*abar2[0]*favg[10]+0.1767766952966369*abar2[1]*favg[6]+0.1767766952966369*abar2[0]*favg[3]; 
  incr[4] = 0.1530931089239486*B0[0]*favg[10]*dv2+0.08838834764831845*B0[1]*favg[6]*dv2+0.08838834764831845*B0[0]*favg[3]*dv2-0.1530931089239486*B1[0]*favg[9]*dv1-0.08838834764831845*B1[1]*favg[5]*dv1-0.08838834764831845*B1[0]*favg[2]*dv1+0.9682458365518543*fjump[14]*amax+0.75*fjump[4]*amax+0.4330127018922193*fjump[0]*amax-0.6846531968814578*abar2[0]*favg[14]-0.3061862178478973*abar2[2]*favg[11]-0.5303300858899107*abar2[1]*favg[8]-0.5303300858899107*abar2[0]*favg[4]-0.3061862178478973*abar2[1]*favg[1]-0.3061862178478973*abar2[0]*favg[0]; 
  incr[5] = (-0.05103103630798287*B0[1]*favg[7]*dv2)+0.1141088661469096*B1[1]*favg[14]*dv1+0.04564354645876385*B1[1]*favg[12]*dv1+0.04564354645876385*B1[1]*favg[11]*dv1+0.0790569415042095*B1[2]*favg[8]*dv1+0.08838834764831845*B1[0]*favg[8]*dv1+0.08838834764831845*B1[1]*favg[4]*dv1+0.04564354645876385*favg[1]*B1[2]*dv1+0.05103103630798287*favg[0]*B1[1]*dv1+0.05103103630798287*B1[0]*favg[1]*dv1-0.25*fjump[5]*amax+0.3061862178478973*abar2[1]*favg[9]+0.158113883008419*abar2[2]*favg[5]+0.1767766952966369*abar2[0]*favg[5]+0.1767766952966369*abar2[1]*favg[2]; 
  incr[6] = (-0.1141088661469096*B0[1]*favg[14]*dv2)-0.04564354645876385*B0[1]*favg[13]*dv2-0.04564354645876385*B0[1]*favg[11]*dv2-0.0790569415042095*B0[2]*favg[8]*dv2-0.08838834764831845*B0[0]*favg[8]*dv2-0.08838834764831845*B0[1]*favg[4]*dv2-0.04564354645876385*favg[1]*B0[2]*dv2-0.05103103630798287*favg[0]*B0[1]*dv2-0.05103103630798287*B0[0]*favg[1]*dv2+0.05103103630798287*B1[1]*favg[7]*dv1-0.25*fjump[6]*amax+0.3061862178478973*abar2[1]*favg[10]+0.158113883008419*abar2[2]*favg[6]+0.1767766952966369*abar2[0]*favg[6]+0.1767766952966369*abar2[1]*favg[3]; 
  incr[7] = (-0.08838834764831845*B0[0]*favg[9]*dv2)-0.05103103630798287*B0[1]*favg[5]*dv2-0.05103103630798287*B0[0]*favg[2]*dv2+0.08838834764831845*B1[0]*favg[10]*dv1+0.05103103630798287*B1[1]*favg[6]*dv1+0.05103103630798287*B1[0]*favg[3]*dv1-0.25*fjump[7]*amax+0.1767766952966369*abar2[0]*favg[7]; 
  incr[8] = 0.1530931089239486*B0[1]*favg[10]*dv2+0.0790569415042095*B0[2]*favg[6]*dv2+0.08838834764831845*B0[0]*favg[6]*dv2+0.08838834764831845*B0[1]*favg[3]*dv2-0.1530931089239486*B1[1]*favg[9]*dv1-0.0790569415042095*B1[2]*favg[5]*dv1-0.08838834764831845*B1[0]*favg[5]*dv1-0.08838834764831845*B1[1]*favg[2]*dv1+0.75*fjump[8]*amax+0.4330127018922193*fjump[1]*amax-0.6846531968814579*abar2[1]*favg[14]-0.2738612787525831*abar2[1]*favg[11]-0.4743416490252572*abar2[2]*favg[8]-0.5303300858899107*abar2[0]*favg[8]-0.5303300858899107*abar2[1]*favg[4]-0.2738612787525831*favg[1]*abar2[2]-0.3061862178478973*abar2[0]*favg[1]-0.3061862178478973*favg[0]*abar2[1]; 
  incr[9] = 0.08838834764831845*B0[0]*favg[7]*dv2-0.1976423537605238*B1[0]*favg[14]*dv1-0.0790569415042095*B1[0]*favg[12]*dv1-0.08838834764831845*B1[2]*favg[11]*dv1-0.1530931089239486*B1[1]*favg[8]*dv1-0.1530931089239486*B1[0]*favg[4]*dv1-0.08838834764831845*favg[1]*B1[1]*dv1-0.08838834764831845*favg[0]*B1[0]*dv1+0.75*fjump[9]*amax+0.4330127018922193*fjump[2]*amax-0.5303300858899107*abar2[0]*favg[9]-0.3061862178478973*abar2[1]*favg[5]-0.3061862178478973*abar2[0]*favg[2]; 
  incr[10] = 0.1976423537605238*B0[0]*favg[14]*dv2+0.0790569415042095*B0[0]*favg[13]*dv2+0.08838834764831845*B0[2]*favg[11]*dv2+0.1530931089239486*B0[1]*favg[8]*dv2+0.1530931089239486*B0[0]*favg[4]*dv2+0.08838834764831845*favg[1]*B0[1]*dv2+0.08838834764831845*favg[0]*B0[0]*dv2-0.08838834764831845*B1[0]*favg[7]*dv1+0.75*fjump[10]*amax+0.4330127018922193*fjump[3]*amax-0.5303300858899107*abar2[0]*favg[10]-0.3061862178478973*abar2[1]*favg[6]-0.3061862178478973*abar2[0]*favg[3]; 
  incr[11] = (-0.08838834764831845*B0[2]*favg[10]*dv2)-0.04564354645876386*B0[1]*favg[6]*dv2-0.05103103630798288*B0[2]*favg[3]*dv2+0.08838834764831845*B1[2]*favg[9]*dv1+0.04564354645876386*B1[1]*favg[5]*dv1+0.05103103630798288*favg[2]*B1[2]*dv1-0.25*fjump[11]*amax+0.3952847075210476*abar2[2]*favg[14]+0.1129384878631565*abar2[2]*favg[11]+0.1767766952966369*abar2[0]*favg[11]+0.2738612787525832*abar2[1]*favg[8]+0.3061862178478973*abar2[2]*favg[4]+0.1767766952966369*favg[0]*abar2[2]+0.158113883008419*abar2[1]*favg[1]; 
  incr[12] = 0.0790569415042095*B1[0]*favg[9]*dv1+0.04564354645876384*B1[1]*favg[5]*dv1+0.04564354645876384*B1[0]*favg[2]*dv1-0.25*fjump[12]*amax+0.1767766952966369*abar2[0]*favg[12]; 
  incr[13] = (-0.0790569415042095*B0[0]*favg[10]*dv2)-0.04564354645876384*B0[1]*favg[6]*dv2-0.04564354645876384*B0[0]*favg[3]*dv2-0.25*fjump[13]*amax+0.1767766952966369*abar2[0]*favg[13]; 
  incr[14] = (-0.1976423537605237*B0[0]*favg[10]*dv2)-0.1141088661469096*B0[1]*favg[6]*dv2-0.1141088661469096*B0[0]*favg[3]*dv2+0.1976423537605237*B1[0]*favg[9]*dv1+0.1141088661469096*B1[1]*favg[5]*dv1+0.1141088661469096*B1[0]*favg[2]*dv1-1.25*fjump[14]*amax-0.9682458365518543*fjump[4]*amax-0.5590169943749475*fjump[0]*amax+0.8838834764831847*abar2[0]*favg[14]+0.3952847075210476*abar2[2]*favg[11]+0.6846531968814578*abar2[1]*favg[8]+0.6846531968814578*abar2[0]*favg[4]+0.3952847075210476*abar2[1]*favg[1]+0.3952847075210476*abar2[0]*favg[0]; 

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

  outl[0] += -1.0*incr[0]*dv12; 
  outl[1] += -1.0*incr[1]*dv12; 
  outl[2] += -1.0*incr[2]*dv12; 
  outl[3] += -1.0*incr[3]*dv12; 
  outl[4] += incr[4]*dv12; 
  outl[5] += -1.0*incr[5]*dv12; 
  outl[6] += -1.0*incr[6]*dv12; 
  outl[7] += -1.0*incr[7]*dv12; 
  outl[8] += incr[8]*dv12; 
  outl[9] += incr[9]*dv12; 
  outl[10] += incr[10]*dv12; 
  outl[11] += -1.0*incr[11]*dv12; 
  outl[12] += -1.0*incr[12]*dv12; 
  outl[13] += -1.0*incr[13]*dv12; 
  outl[14] += -1.0*incr[14]*dv12; 
return std::abs(amid); 
} 
